function build_residual!(residual::Matrix{Float64}, u::Matrix{Float64}, t::Float64, BChandler::Dict, dg::DGStd, param::parameters)
    # For a linear mesh, we can simplify the computation
    if dg.mesh isa LMesh
        # We compute the projected reference flux
        flux = compute_physflux(block_matmul(dg.refelem.chiq, u, dg.mesh.Nel), param)
        flux_to_ref!(flux, dg)
        flux = Tuple(block_matmul(dg.refelem.Ph, flux[dir], dg.mesh.Nel) for dir in 1:dg.dim)

        # Now, we interpolate the projected flux to the faces
        fluxface = Tuple(block_matmul(dg.refelem.chif, flux[dir], dg.mesh.Nel) for dir in 1:dg.dim)

        # Finally, we evaluate numflux
        un = block_matmul(dg.refelem.chif, u, dg.mesh.Nel)
        up = dg.FtoF * un + evaluate_BC(BChandler, dg, t)
        numflux = compute_numflux(un, up, dg.nphys, param)
        flux_to_ref!(numflux, dg)

        # Assemble complete residual (volume and face)
        residual .= 0.0
        for dir in 1:dg.dim
            @views residual .= residual .- block_matmul(dg.refelem.Dh[dir], flux[dir], dg.mesh.Nel) .- block_matmul((dg.refelem.LIFT[dir]), numflux[dir] .- fluxface[dir], dg.mesh.Nel)
        end

        # Scale by Jacobian
        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            @views residual[index,:] .= residual[index,:] ./ dg.mesh.detJ[ielem]
        end

        # Add source (if applicable)
        if !(isnothing(param.sourcename))
            residual .= residual .+ block_matmul(dg.refelem.Ph, compute_source(dg, param, dg.qpts, t), dg.mesh.Nel)
        end

        return residual
    end
end

function build_residual!(residual::Matrix{Float64}, u::Matrix{Float64}, t::Float64, BChandler::Dict, dg::DGFluxDiff, param::parameters)
    if dg.mesh isa LMesh
        # We start by computing the entropy-projected solution (volume and face)
        v = compute_evar(block_matmul(dg.refelem.chiq, u, dg.mesh.Nel), param) # Compute entropy variables at vol quadrature points
        v = block_matmul(dg.refelem.Ph, v, dg.mesh.Nel) # Project entropy variables
        vq = block_matmul(dg.refelem.chiq, v, dg.mesh.Nel) # evaluate at vol quadrature
        vf = block_matmul(dg.refelem.chif, v, dg.mesh.Nel)

        uq = compute_cvar(vq, param)
        un = compute_cvar(vf, param)

        up = dg.FtoF * un + evaluate_BC(BChandler, dg, t)

        vf = vq = v = nothing
        residual .= 0.0

        # Volume terms
        Npts = dg.refelem.Nfnodes*dg.refelem.Nfaces + dg.refelem.Nqnodes # number of 2-pt flux pts
        F = Tuple(Array{Float64}(undef, Npts, Npts, dg.Nstates) for dir in 1:dg.dim)

        @inbounds for dir in 1:param.dim # Allocate diagonals with zeros (don't need them since Hadamard prod with Skew-symmetric)
            for j in 1:Npts
                F[dir][j,j,:] .= 0.0
            end
        end

        for ielem = 1:dg.mesh.Nel
            indexv = 1+dg.refelem.Nqnodes*(ielem-1):dg.refelem.Nqnodes*ielem
            indexf = 1+dg.refelem.Nfnodes*dg.refelem.Nfaces*(ielem-1):dg.refelem.Nfnodes*dg.refelem.Nfaces*ielem
            indexb = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            
            uv = @view uq[indexv,:]
            uf = @view un[indexf,:]

            F = compute_two_pt_flux!(F, uv, uf, param) # in place computation two-point flux matrix
            two_pt_flux_to_ref!(F, ielem, dg)

            for dir in 1:dg.dim
                @views residual[indexb,:] .= residual[indexb,:] .- dg.refelem.MVF * reshape(sum(dg.refelem.SS[dir] .* F[dir], dims=2), (Npts, dg.Nstates))
            end
        end

        # Surface term
        numflux = compute_numflux(un, up, dg.nphys, param)
        flux_to_ref!(numflux, dg)
        for dir in 1:dg.dim
            @views residual .= residual .- block_matmul((dg.refelem.LIFT[dir]), numflux[dir], dg.mesh.Nel)
        end

        # Scale by Jacobian
        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            @views residual[index,:] .= residual[index,:] ./ dg.mesh.detJ[ielem]
        end

        # Add source (if applicable)
        if !(isnothing(param.sourcename))
            residual .= residual .+ block_matmul(dg.refelem.Ph, compute_source(dg, param, dg.qpts, t), dg.mesh.Nel)
        end

        return residual
    end
end

function build_residual!(residual::Matrix{Float64}, u::Matrix{Float64}, t::Float64, BChandler::Dict, dg::DGArtVisc, param::parameters)
end

function build_residual!(residual::Matrix{Float64}, u::Matrix{Float64}, t::Float64, BChandler::Dict, dg::DGEntFilt, param::parameters)
end

function block_matmul(block::AbstractArray, myvec::AbstractArray, N::Integer) # Block * v
    Nrowv = size(myvec,1)
    Ncolv = size(myvec,2)
    Nrowb = size(block,1)
    Ncolb = size(block,2)

    return (reshape(block * reshape(myvec, (Ncolb,div(Nrowv * Ncolv,Ncolb))), (N*Nrowb, Ncolv)))
end

function block_matmul!(out::AbstractArray, block::AbstractArray, myvec::AbstractArray, N::Integer)
    Nrowv = size(myvec,1)
    Ncolv = size(myvec,2)
    Nrowb = size(block,1)
    Ncolb = size(block,2)

    mul!(out, block, reshape(myvec, (Ncolb,div(Nrowv * Ncolv,Ncolb))))
    out = reshape(out, (N*Nrowb, Ncolv))

    return nothing
end

function block_matmul_add!(out::AbstractArray, block::AbstractArray, myvec::AbstractArray, N::Integer)
    Nrowv = size(myvec,1)
    Ncolv = size(myvec,2)
    Nrowb = size(block,1)
    Ncolb = size(block,2)

    mul!(reshape(out, (Nrowb, div(Nrowv * Ncolv,Ncolb))), block, reshape(myvec, (Ncolb,div(Nrowv * Ncolv,Ncolb))), 1.0, 1.0)
    out = reshape(out, (N*Nrowb, Ncolv))

    return nothing
end

# OPTIMIZE THIS FUNCTION LATER. MAYBE IN-PLACE REPLACEMENT IS NOT THE BEST OPTION (I was only thinking about)...
function flux_to_ref!(f::Tuple, dg::DG) # this is a global operation
    if dg.dim == 1
        return nothing
    elseif dg.dim == 2
        if dg.mesh isa LMesh
            ftemp = zeros((dg.refelem.Nbnodes, dg.Nstates))
            for ielem in 1:dg.mesh.Nel
                index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
                @views ftemp .= dg.mesh.CT[ielem][1,1] .* f[1][index,:] .+ dg.mesh.CT[ielem][1,2] .* f[2][index,:]
                @views f[2][index,:] .= dg.mesh.CT[ielem][2,1] .* f[1][index,:] .+ dg.mesh.CT[ielem][2,2] .* f[2][index,:]
                @views f[1][index,:] .= ftemp
            end
            return nothing
        end
    end
end

# OPTIMIZE THIS FUNCTION LATER. MAYBE IN-PLACE REPLACEMENT IS NOT THE BEST OPTION...
function two_pt_flux_to_ref!(F::Tuple, iele::Integer, dg::DG) # this is a local operation
    if dg.dim == 1
        return nothing
    elseif dg.dim == 2
        if dg.mesh isa LMesh
            Ftemp = zeros(size(F))
            @views Ftemp .= dg.mesh.CT[ielem][1,1] .* F[1] .+ dg.mesh.CT[ielem][1,2] .* F[2]
            @views F[2] .= dg.mesh.CT[ielem][2,1] .* F[1] .+ dg.mesh.CT[ielem][2,2] .* F[2]
            @views F[1] .= Ftemp
            return nothing
        end
    end
end
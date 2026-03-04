#####################################################################
# Assemble residual for different DG flavours
#####################################################################

function build_residual!(residual::Matrix{Float64}, u::Matrix{Float64}, t::Float64, BChandler::Dict, dg::DGStd, param::parameters)
    # For a linear mesh, we can simplify the computation
    if dg.mesh isa LMesh
        # We compute the projected reference flux
        flux = compute_physflux(block_matmul(dg.refelem.chiq, u, dg.mesh.Nel), param)
        flux_to_ref!(flux, dg.refelem.Nqnodes, dg)
        flux = Tuple(block_matmul(dg.refelem.Ph, flux[dir], dg.mesh.Nel) for dir in 1:dg.dim)

        # Now, we interpolate the projected flux to the faces
        fluxface = Tuple(block_matmul(dg.refelem.chif, flux[dir], dg.mesh.Nel) for dir in 1:dg.dim)

        # Finally, we evaluate numflux
        un = block_matmul(dg.refelem.chif, u, dg.mesh.Nel)
        up = dg.FtoF * un + evaluate_BC(BChandler, dg, t)
        numflux = compute_numflux(un, up, dg.nphys, param)
        flux_to_ref!(numflux, dg.refelem.Nfnodes, dg)

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
        flux_to_ref!(numflux, dg.refelem.Nqnodes, dg)
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
    if dg.mesh isa LMesh
        # We start by computing the projected entropy variables and we evaluate face quantities
        uq = block_matmul(dg.refelem.chiq, u, dg.mesh.Nel)
        v = compute_evar(uq, param)
        v = block_matmul(dg.refelem.Ph, v, dg.mesh.Nel)

        vn = block_matmul(dg.refelem.chif, v, dg.mesh.Nel)
        un = compute_cvar(vn, param)
        up = dg.FtoF * un + evaluate_BC(BChandler, dg, t)
        vp = compute_evar(up, param) # FIX FOR BOUNDARY CONDITIONS!

        # FIRST AUXILIARY PROBLEM
        theta = AV_auxiliary1(v, vn, vp, dg)

        # We compute projected scaled entropy gradient
        thetaq = Tuple(block_matmul(dg.refelem.chiq, theta[dir], dg.mesh.Nel) for dir in 1:dg.dim)
        thetaK = Tuple(Matrix{Float64}(undef, dg.mesh.Nel*dg.refelem.Nqnodes, dg.Nstates) for dir in 1:dg.dim)

        for iqnode in 1:dg.mesh.Nel*dg.refelem.Nqnodes
            K = compute_cvar_Hessian(uq[iqnode,:], param)
            for dir in 1:dg.dim
                @views thetaK[dir][iqnode,:] .= K * thetaq[dir][iqnode,:]
            end
        end
        thetaK = Tuple(block_matmul(dg.refelem.Ph, thetaK[dir], dg.mesh.Nel) for dir in 1:dg.dim)

        # We compute the projected reference flux (we'll need it a few times)
        flux = compute_physflux(uq, param)
        flux_to_ref!(flux, dg.refelem.Nqnodes, dg)
        flux = Tuple(block_matmul(dg.refelem.Ph, flux[dir], dg.mesh.Nel) for dir in 1:dg.dim)

        # Artificial viscosity
        sigma = Tuple(Matrix{Float64}(undef, dg.DOF, dg.Nstates) for dir in 1:dg.dim)
        for ielem = 1:dg.mesh.Nel
            indexf = 1+dg.refelem.Nfnodes*dg.refelem.Nfaces*(ielem-1):dg.refelem.Nfnodes*dg.refelem.Nfaces*ielem
            indexb = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem

            # Elemental entropy deficit
            @views psi = compute_cvar_potential(un[indexf,:], param)
            flux_to_ref!(psi, dg.refelem.Nfnodes, dg)

            delta = 0.0
            for dir in 1:dg.dim
                for istate in 1:dg.Nstates
                    @views delta -= flux[dir][indexb,istate]' * dg.refelem.Qh[dir] * v[indexb,istate]
                end
                @views delta += dot(dg.refelem.bh[dir], psi[dir])
            end

            # Entropy denominator
            den = 0.0
            for dir in 1:dg.dim, istate in 1:dg.Nstates
                @views den += sum(thetaK[dir][indexb,istate]' * dg.refelem.M * theta[dir][indexb,istate]) # since we are using a block diagonal K as in (Chan 2025)
            end
            den *= dg.mesh.detJ[ielem] # scale by Jacobian

            # Artifical viscosity coefficient
            epsilon = AV_coeff(delta, den, param)

            # We finally build the viscous entropy fluxes
            for dir in 1:dg.dim
                @views sigma[dir][indexb,:] .= epsilon .* thetaK[dir][indexb,:]
            end
        end

        # THIRD AUXILIARY PROBLEM
        sigman = Tuple(block_matmul(dg.refelem.chif, sigma[dir], dg.mesh.Nel) for dir in 1:dg.dim)
        sigmap = Tuple(dg.FtoF * sigman[dir] for dir in 1:dg.dim) # FIX FOR BOUNDARY CONDITIONS
        gvisc = AV_auxiliary2(sigma, sigman, sigmap, dg)

        # PRIMARY PROBLEM
        numflux = compute_numflux(un, up, dg.nphys, param)
        flux_to_ref!(numflux, dg.refelem.Nfnodes, dg)

        residual .= 0.0
        for dir in 1:dg.dim
            @views residual .= residual .+ block_matmul(dg.refelem.MinvQhT[dir], flux[dir], dg.mesh.Nel) .- block_matmul((dg.refelem.LIFT[dir]), numflux[dir], dg.mesh.Nel) # weak DG
        end

        # Scale by Jacobian
        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            @views residual[index,:] .= residual[index,:] ./ dg.mesh.detJ[ielem]
        end

        # Viscous correction
        residual .= residual .+ gvisc
        
        # Add source (if applicable)
        if !(isnothing(param.sourcename))
            residual .= residual .+ block_matmul(dg.refelem.Ph, compute_source(dg, param, dg.qpts, t), dg.mesh.Nel)
        end

        return residual
    end
end

function build_residual!(residual::Matrix{Float64}, u::Matrix{Float64}, t::Float64, BChandler::Dict, dg::DGEntFilt, param::parameters)
end

#####################################################################
# Other helper functions
#####################################################################

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
function flux_to_ref!(f::Tuple, Nnodes::Integer, dg::DG) # this is a global operation
    if dg.dim == 1
        return nothing
    elseif dg.dim == 2
        if dg.mesh isa LMesh
            ftemp = zeros((Nnodes, dg.Nstates))
            for ielem in 1:dg.mesh.Nel
                index = 1+Nnodes*(ielem-1):Nnodes*ielem
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

function grad_to_ref(u::AbstractArray, Nnodes::Integer, dir::Integer, dg::DG) # this is a global operation
    if dg.dim == 1
        return (u,)
    elseif dg.dim == 2
        if dg.mesh isa LMesh
            f = (Matrix{Float64}(undef, Nnodes, dg.Nstates), Matrix{Float64}(undef, Nnodes, dg.Nstates))
            for ielem in 1:dg.mesh.Nel
                index = 1+Nnodes*(ielem-1):dg.Nnodes*ielem
                f[1][index,:] .= dg.mesh.CT[1,dir] .* u
                f[2][index,:] .= dg.mesh.CT[2,dir] .* u
            end
            return f
        end
    end
end

#####################################################################
# Auxiliary functions for AV
#####################################################################

# MUST USE PROJECTED ENTROPY VARIABLES AS INPUT
function AV_auxiliary1(v::AbstractMatrix, vn::AbstractMatrix, vp::AbstractMatrix, dg::DGArtVisc)
    if dg.mesh isa LMesh
        theta = Tuple(zeros(dg.DOF, dg.Nstates) for dir in 1:dg.dim)
        vjump = 0.5 .* (vp .- vn) # We use a penalty term as in (Chan, 2025)

        for dir1 in 1:dg.dim
            # Map entropy quantities to ref elem
            v = grad_to_ref(v, dg.refelem.Nbnodes, dir1, dg)
            vjump = grad_to_ref(vjump, dg.refelem.Nbnodes, dir1, dg)

            for dir2 in 1:dg.dim
            @views theta[dir1] .= theta[dir1] .+ block_matmul(dg.refelem.Dh[dir2], v[dir2], dg.mesh.Nel) .+ 0.5 .* block_matmul((dg.refelem.LIFT[dir2]), vjump[dir2], dg.mesh.Nel)
            end

            # We finally scale by Jacobian
            for ielem = 1:dg.mesh.Nel
                index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
                @views theta[dir1][index,:] .= theta[dir1][index,:] ./ dg.mesh.detJ[ielem]
            end
        end

        return theta
    end
end

function AV_auxiliary2(sigma::Tuple, sigman::Tuple, sigmap::Tuple, dg::DGArtVisc)
    if dg.mesh isa LMesh
        gvisc = zeros(dg.DOF, dg.Nstates)
        sigmanum = 0.5 .* (sigman .+ sigmap) # We use central viscous fluxes as in (Chan, 2025)
        flux_to_ref!(sigma, dg.refelem.Nbnodes, dg)
        flux_to_ref!(sigmanum, dg.refelem.Nfnodes, dg)

        for dir in 1:dg.dim
            @views gvisc .= gvisc .- block_matmul(dg.refelem.MinvQhT[dir], sigma[dir], dg.mesh.Nel) .+ block_matmul((dg.refelem.LIFT[dir]), sigmanum[dir], dg.mesh.Nel)
        end

        # We finally scale by Jacobian
        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            @views gvisc[index,:] .= gvisc[index,:] ./ dg.mesh.detJ[ielem]
        end

        return gvisc
    end
end
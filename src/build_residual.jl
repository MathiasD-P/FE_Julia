function build_residual(u, t, BChandler::Dict, solver::SolverStd, dg::DGStd, param::parameters)
    # For a linear mesh, we can simplify the computation
    if dg.mesh isa LMesh
        # We compute the projected reference flux
        uq = block_matmul!(solver.uq, dg.refelem.chiq, u, dg.mesh.Nel)
        fluxq = compute_physflux!(solver.fluxq, uq, param)
        flux_to_ref!(fluxq, dg)

        for dir in 1:dg.dim
            @views flux = block_matmul!(solver.flux[dir], dg.refelem.Ph, solver.fluxq[dir], dg.mesh.Nel)
        end

        # Now, we interpolate the projected flux to the faces
        for dir in 1:dg.dim
            @views fluxface = block_matmul!(solver.fluxface[dir], dg.refelem.chif, solver.flux[dir], dg.mesh.Nel)
        end

        # Finally, we evaluate numflux
        un = block_matmul!(solver.un, dg.refelem.chif, u, dg.mesh.Nel)
        up = block_matmul!(solver.up, dg.FtoF, un, dg.mesh.Nel)
        up .= up .+ evaluate_BC(BChandler, dg, t)

        
        numflux = compute_numflux(solver.numflux, un, up, dg.nphys, param)
        flux_to_ref!(numflux, dg)

        # Assemble complete residual (volume and face)
        solver.residual .= 0.0
        residual = solver.residual
        for dir in 1:dg.dim
            @views residual .= residual .- block_matmul!(solver.resbuffer1, dg.refelem.Dh[dir], solver.flux[dir], dg.mesh.Nel) .- block_matmul!(solver.resbuffer2,dg.refelem.LIFT[dir], numflux[dir] .- solver.fluxface[dir], dg.mesh.Nel)
        end

        # Scale by Jacobian
        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            @views residual[index,:] .= residual[index,:] ./ dg.mesh.J[ielem]
        end

        return residual
    end
end

function build_residual(u, t, BChandler::Dict, solver, dg::DGFluxDiff, param::parameters)
    if dg.mesh isa LMesh
        # We start by computing the entropy-projected solution (volume and face)
        v = eval_evar(block_matmul(dg.refelem.chiq, u, dg.mesh.Nel), param) # Compute entropy variables at vol quadrature points
        v = block_matmul(dg.refelem.Ph, v, dg.mesh.Nel) # Project entropy variables
        vq = block_matmul(dg.refelem.chiq, v, dg.mesh.Nel) # evaluate at vol quadrature
        vf = block_matmul(dg.refelem.chif, v, dg.mesh.Nel)

        uq = eval_cvar(vq, param)
        un = eval_cvar(vf, param)
        up = dg.FtoF * un + evaluate_BC(BChandler, dg, t)

        residual = zeros(size(u))

        # Volume terms
        for ielem = 1:dg.mesh.Nel
            Npts = dg.refelem.Nfnodes*dg.refelem.Nfaces + dg.refelem.Nqnodes
            indexv = 1+dg.refelem.Nqnodes*(ielem-1):dg.refelem.Nqnodes*ielem
            indexf = 1+dg.refelem.Nfnodes*dg.refelem.Nfaces*(ielem-1):dg.refelem.Nfnodes*dg.refelem.Nfaces*ielem
            indexb = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            
            uv = @view uq[indexv,:]
            uf = @view un[indexf,:]

            F = two_pt_flux(uv, uf, param)
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
            @views residual[index,:] .= residual[index,:] ./ dg.mesh.J[ielem]
        end

        return residual
    end
end

function build_residual(u, t, BChandler::Dict, dg::DGArtVisc, param::parameters)
end

function build_residual(u, t, BChandler::Dict, dg::DGEntFilt, param::parameters)
end

function block_matmul(block, myvec, N::Integer) # Block * v
    Nrowv = size(myvec,1)
    Ncolv = size(myvec,2)
    Nrowb = size(block,1)
    Ncolb = size(block,2)

    return (reshape(block * reshape(myvec, (Ncolb,div(Nrowv * Ncolv,Ncolb))), (N*Nrowb, Ncolv)))
end

function block_matmul!(out, block, myvec, N::Integer)
    Nrowv = size(myvec,1)
    Ncolv = size(myvec,2)
    Nrowb = size(block,1)
    Ncolb = size(block,2)
    Ntemp = div(Nrowv * Ncolv,Ncolb)

    reshapedview = reshape(out, (Nrowb,Ntemp))
    mul!(reshapedview, block, reshape(myvec, (Ncolb,Ntemp)))

    return out
end

# OPTIMIZE THIS FUNCTION LATER. MAYBE IN-PLACE REPLACEMENT IS NOT THE BEST OPTION (I was only thinking about)...
function flux_to_ref!(f, dg::DG) # this is a global operation
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
function two_pt_flux_to_ref!(F, ielem, dg::DG) # this is a local operation
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
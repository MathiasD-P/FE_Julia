function build_residual(u, t, BChandler::Dict, dg::DGStd, param::parameters)
    # For a linear mesh, we can simplify the computation
    if dg.mesh isa LMesh
        # We compute the projected reference flux
        flux = compute_physflux(block_matmul(dg.refelem.chiq, u, dg.mesh.Nel), param)
        flux_to_ref!(flux, dg)
        for dir in 1:dg.dim
            @views flux[dir] = block_matmul(dg.refelem.Ph, flux[dir], dg.mesh.Nel)
        end

        # Now, we interpolate the projected flux to the faces
        fluxface = Vector{Matrix{Float64}}(undef, dg.dim)
        for dir in 1:dg.dim
            @views fluxface[dir] = block_matmul(dg.refelem.chif, flux[dir], dg.mesh.Nel)
        end

        # Finally, we evaluate numflux
        un = block_matmul(dg.refelem.chif, u, dg.mesh.Nel)
        up = dg.FtoF * un .+ evaluate_BC(BChandler, dg, t)
        numflux = compute_numflux(un, up, dg.nphys, param)
        flux_to_ref!(numflux, dg)

        # Assemble complete residual
        residual = zeros(size(u))
        for dir in 1:dg.dim
            residual .= residual .- block_matmul(dg.refelem.Dh[dir], flux[dir], dg.mesh.Nel) .- block_matmul((dg.refelem.LIFT[dir]), numflux[dir] .- fluxface[dir], dg.mesh.Nel)
        end

        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
            @views residual[index,:] .= residual[index,:] ./ dg.mesh.J[ielem]
        end

        return residual
    end
end

function build_residual(u, t, BChandler::Dict, dg::DGFluxDiff, param::parameters)
end

function build_residual(u, t, BChandler::Dict, dg::DGArtVisc, param::parameters)
end

function build_residual(u, t, BChandler::Dict, dg::DGEntFilt, param::parameters)
end

function block_matmul(block, myvec, N::Integer) # Block * v
    # prod = Matrix{Float64}(undef, size(block,1) * N, size(vec,2))
    # @inbounds @simd for i = 1:N
    #     prod[size(block,1)*(i-1)+1:size(block,1)*i,:] = block * vec[size(block,2)*(i-1)+1:size(block,2)*i,:]
    # end

    # return prod


    Nrowv, Ncolv = size(myvec)
    Nrowb, Ncolb = size(block)

    return (reshape(block * reshape(myvec, (Ncolb,div(Nrowv * Ncolv,Ncolb))), (N*Nrowb, Ncolv)))
end

function block_matmul!(block, myvec, out, N::Integer)
    Nrowv, Ncolv = size(myvec)
    Nrowb, Ncolb = size(block)

    mul!(out, block, reshape(myvec, (Ncolb,div(Nrowv * Ncolv,Ncolb))))
    out = reshape(out, (N*Nrowb, Ncolv))

    return nothing
end

function flux_to_ref!(f, dg::DGStd)
    if dg.dim == 1
        return nothing
    elseif dg.dim == 2
        if dg.mesh isa LMesh
            for i in 1:dg.dim
                for ielem in 1:dg.mesh.Nel
                    index = 1+dg.refelem.Nbnodes*(ielem-1):dg.refelem.Nbnodes*ielem
                    @views f[i][index,:] .= dg.mesh.CT[ielem][i,1] .* f[1][index,:] .+ dg.mesh.CT[ielem][i,2] .* f[2][index,:]
                end
            end
            return nothing
        end
    end
end
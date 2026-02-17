#####################################################################
# Everything we want to do with the solution once we have it
#####################################################################

function compute_L2error(u, t, enodes::AbstractNodes, dg::DG, param::parameters)
    if isa(dg.mesh, LMesh)
        if typeof(enodes) != typeof(dg.refelem.bnodestype)
            error("Node type used for the computation of the L2 error does not agree with the basis nodes!")
        end

        Nenodes = numnodes(enodes)

        chie, we, epts = extract_volume_quadrature(dg.refelem.bnodestype, enodes)
        epts = reduce(vcat, mapping(dg.mesh, epts, ielem) for ielem in 1:dg.mesh.Nel)

        if param.pdetype == "LinAdv"
            if param.BCname == "periodic"
                if param.domain == "unit_interval_linear"
                    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, mod.((epts .- param.a .* t .+ 0.5),1.0) .- 0.5)).^2
                end
            end
        
        elseif param.pdetype == "EulerPerfGas"
            if param.BCname == "periodic"
                if param.ICname == "IsentropicDensityWave"
                    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, epts .- 0.1 * t)).^2
                
                elseif param.ICname == "GassnerEuler"
                    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, epts .- 2 * t)).^2
                end
            end
        end

        for ielem = 1:dg.mesh.Nel
            index = 1+Nenodes*(ielem-1):Nenodes*ielem
            @views error2[index,:] .= error2[index,:] .* dg.mesh.J[ielem]
        end

        return sum(sqrt.(sum(block_matmul(Diagonal(we), error2, dg.mesh.Nel), dims=1)))
    end
end
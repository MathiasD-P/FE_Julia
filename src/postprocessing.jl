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
        
        elseif param.pdetype == "Burgers"
            if param.BCname == "periodic"
                if param.ICname == "GassnerBurgers"
                    c = 2.2
                    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, epts .- c * t)).^2

                elseif param.ICname == "sin_1state"
                    f = ones(size(epts))
                    utrue = ones(size(epts))
                    i = 0
                    while maximum(abs.(f)) > 1e-13
                        fsin!(param.k[1], param.phi[1], param.av, epts, t, utrue, f)
                        Newton_step_sin!(param.k[1], param.phi[1], epts, t, utrue, f)
                        i += 1

                        if i > 2000
                            error("Newton iteration did not converge!")
                        end
                    end

                    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- utrue).^2
                end

            elseif param.BCname == "unit_shock"
                if param.ICname == "tanh_1state"
                    f = ones(size(epts))
                    utrue = ones(size(epts))
                    i = 0
                    while maximum(abs.(f)) > 1e-13
                        ftanh!(param.k[1], param.phi[1], epts, t, utrue, f)
                        Newton_step_tanh!(param.k[1], param.phi[1], epts, t, utrue, f)
                        i += 1

                        if i > 2000
                            error("Newton iteration did not converge!")
                        end
                    end

                    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- utrue).^2
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

# HELPER FUNCTIONS FOR IMPLICIT BURGERS SOLVE
function fsin!(k, phi, av, x, t, u, dest)
    dest .= sin.(2*pi*k .* (x .- t .* u .- phi)) .- u .+ av
end

function Newton_step_sin!(k, phi, x, t, u, f)
    u .= u .- f ./ (-2*pi*t*k .* cos.(2*pi*k .* (x .- t .* u .- phi)) .- 1.0)
end

function ftanh!(k, phi, x, t, u, dest)
    dest .= tanh.(k .* (x .- t .* u .- phi)) .- u
end

function Newton_step_tanh!(k, phi, x, t, u, f)
    u .= u .- f ./ (-k*t .* sech.(k .* (x .- t .* u .- phi)).^2 .- 1.0)
end
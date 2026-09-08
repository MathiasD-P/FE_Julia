#####################################################################
# Everything we want to do with the solution once we have it
#####################################################################

function local_L2error(u, t, epts, chie, dg::DG, PDE::LinAdv, ic::InitialCondition, BChandler::Dict, source::Nothing)
    if !(isperiodic(BChandler) && isperiodicunitinterval(dg.mesh))
        error("Error for linear advection can only be computed for periodic BCs on the unit interval!")
    end

    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(ic, dg, PDE, mod.((epts .- PDE.a .* t .+ 0.5),1.0) .- 0.5)).^2
    return error2
end

function local_L2error(u, t, epts, chie, dg::DG, PDE::Burgers, ic::GassnerBurgers, BChandler::Dict, source::GassnerBurgersSource)
    if !(isperiodic(BChandler) && isperiodicunitinterval(dg.mesh))
        error("Error for Burgers can only be computed for periodic BCs on the unit interval!")
    end

    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(ic, dg, PDE, epts .- source.c * t)).^2
    return error2
end

function local_L2error(u, t, epts, chie, dg::DG, PDE::Burgers, ic::Sin1State, BChandler::Dict, source::Nothing)
    if !(isperiodic(BChandler) && isperiodicunitinterval(dg.mesh))
        error("Error for Burgers can only be computed for periodic BCs on the unit interval!")
    end

    if t > 1 / ic.k - 0.0001
        error("The time selected is too close to shock formation!")
    end

    # We use a Newton iteration to compute the exact solution from the implict equation u(x,t) = f(x-u(x,t)t)

    tol = 1e-13 # tolerance for Newton iteration
    MaxIter = 2000 # maximum number of iterations for Newton iteration

    # Helper functions for Newton iteration!
    function fsin!(k, phi, av, x, t, ucurr, dest)
        dest .= sin.(k .* (x .- t .* ucurr .- phi)) .- ucurr .+ av
    end

    function Newton_step_sin!(k, phi, x, t, ucurr, f)
        ucurr .= ucurr .- f ./ (-t*k .* cos.(k .* (x .- t .* ucurr .- phi)) .- 1.0)
    end

    f = ones(size(epts))
    utrue = ones(size(epts))
    i = 0
    while maximum(abs.(f)) > tol
        fsin!(ic.k, ic.phi, ic.av, epts, t, utrue, f)
        Newton_step_sin!(ic.k, ic.phi, epts, t, utrue, f)
        i += 1

        if i > MaxIter
            error("Newton iteration did not converge!")
        end
    end

    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- utrue).^2
    return error2
end

function local_L2error(u, t, epts, chie, dg::DG, PDE::EulerPerfGas, ic::GassnerEuler, BChandler::Dict, source::GassnerEulerSource)
    if !(isperiodic(BChandler) && isperiodicunitinterval(dg.mesh))
        error("Error for Euler can only be computed for periodic BCs on the unit interval!")
    end

    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(ic, dg, PDE, epts .- 2 * t)).^2
    return error2
end

function local_L2error(u, t, epts, chie, dg::DG, PDE::EulerPerfGas, ic::IsentropicDensityWave, BChandler::Dict, source::Nothing)
    if !(isperiodic(BChandler) && isperiodicunitinterval(dg.mesh))
        error("Error for Euler can only be computed for periodic BCs on the unit interval!")
    end

    error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(ic, dg, PDE, epts .- 0.1 * t)).^2
    return error2
end


function compute_L2error(u, t, enodes::AbstractNodes, ic::InitialCondition, dg::DG, physics::PhysProp)
     if isa(dg.mesh, LMesh)
        if typeof(enodes) != typeof(dg.refelem.bnodestype)
            error("Node type used for the computation of the L2 error does not agree with the basis nodes!")
        end

        Nenodes = numnodes(enodes)

        chie, we, epts = extract_volume_quadrature(dg.refelem.bnodestype, enodes)
        epts = reduce(vcat, mapping(dg.mesh, epts, ielem) for ielem in 1:dg.mesh.Nel)

        error2 = local_L2error(u, t, epts, chie, dg, physics.PDE, ic, physics.BChandler, physics.source)

        for ielem = 1:dg.mesh.Nel
            index = 1+Nenodes*(ielem-1):Nenodes*ielem
            @views error2[index,:] .= error2[index,:] .* dg.mesh.J[ielem]
        end

        return sum(sqrt.(sum(block_matmul(Diagonal(we), error2, dg.mesh.Nel), dims=1)))
    end
end

# function compute_L2error(u, t, enodes::AbstractNodes, dg::DG, param::parameters)
#     if isa(dg.mesh, LMesh)
#         if typeof(enodes) != typeof(dg.refelem.bnodestype)
#             error("Node type used for the computation of the L2 error does not agree with the basis nodes!")
#         end

#         Nenodes = numnodes(enodes)

#         chie, we, epts = extract_volume_quadrature(dg.refelem.bnodestype, enodes)
#         epts = reduce(vcat, mapping(dg.mesh, epts, ielem) for ielem in 1:dg.mesh.Nel)

#         if param.pdetype == "LinAdv"
#             if param.BCname == "periodic"
#                 if param.domain == "unit_interval_linear"
#                     error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, mod.((epts .- param.a .* t .+ 0.5),1.0) .- 0.5)).^2
#                 end
#             end
        
#         elseif param.pdetype == "Burgers"
#             if param.BCname == "periodic"
#                 if param.ICname == "GassnerBurgers"
#                     c = 2.2
#                     error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, epts .- c * t)).^2

#                 elseif param.ICname == "sin_1state"
#                     f = ones(size(epts))
#                     utrue = ones(size(epts))
#                     i = 0
#                     while maximum(abs.(f)) > 1e-13
#                         fsin!(param.k, param.phi, param.av, epts, t, utrue, f)
#                         Newton_step_sin!(param.k, param.phi, epts, t, utrue, f)
#                         i += 1

#                         if i > 2000
#                             error("Newton iteration did not converge!")
#                         end
#                     end

#                     error2 = (block_matmul(chie, u, dg.mesh.Nel) .- utrue).^2
#                 end
#             end
        
#         elseif param.pdetype == "EulerPerfGas"
#             if param.BCname == "periodic"
#                 if param.ICname == "IsentropicDensityWave"
#                     error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, epts .- 0.1 * t)).^2
                
#                 elseif param.ICname == "GassnerEuler"
#                     error2 = (block_matmul(chie, u, dg.mesh.Nel) .- initialize_states(dg, param, epts .- 2 * t)).^2
#                 end
#             end
#         end

#         for ielem = 1:dg.mesh.Nel
#             index = 1+Nenodes*(ielem-1):Nenodes*ielem
#             @views error2[index,:] .= error2[index,:] .* dg.mesh.J[ielem]
#         end

#         return sum(sqrt.(sum(block_matmul(Diagonal(we), error2, dg.mesh.Nel), dims=1)))
#     end
# end

# # HELPER FUNCTIONS FOR IMPLICIT BURGERS SOLVE
# function fsin!(k, phi, av, x, t, u, dest)
#     dest .= sin.(2*pi*k .* (x .- t .* u .- phi)) .- u .+ av
# end

# function Newton_step_sin!(k, phi, x, t, u, f)
#     u .= u .- f ./ (-2*pi*t*k .* cos.(2*pi*k .* (x .- t .* u .- phi)) .- 1.0)
# end
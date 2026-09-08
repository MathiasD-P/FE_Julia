# function run(param_file)
#     param = parse_parameters(param_file)

#     if param.OOAtest
#         Nrefinements = param.Nrefinements
#         L2error = zeros((Nrefinements,)) # store DOFS and Nrefinements
#     else
#         Nrefinements = 1
#     end

#     for isolve = 1:Nrefinements
#         if param.OOAtest
#             sol, L2error[isolve,:] = set_up_and_solve(param)
#             refinemesh!(param)
#         else
#             sol = set_up_and_solve(param)
#         end

#         N = length(sol)

#         if N == 2
#             outputnames = ("solution", "time")
#         elseif N == 3
#             outputnames = ("solution", "time", "entropy")
#         end

#         for iname in 1:N
#             if param.save
#             else
#                 println(name[i] * " :")
#                 println(sol[])
#             end
#         end
#     end

#     if param.save
#     else
#         println("L2Error: ")
#         println(L2error)
#     end

# end

function set_up_problem(param::parameters)
    # Initialize Problem Physics
    PDE = parse_parameters_PDE(param)
    ic, source = parse_parameters_initialcondition_and_source(param)

    numflux = parse_parameters_numflux(param)
    tpflux = parse_parameters_tpflux(param)
    artviscmodel = parse_parameters_artvisc(param)
    rescorrmodel = parse_parameters_rescorr(param)


    # initialize mesh and nodes
    mesh = initialize_mesh(param)
    bnodes = make_nodes(param.bnodes)
    qnodes = make_nodes(param.qnodes)
    fnodes = make_nodes(param.refelem, param.fnodes)

    # initialize ref element and DG object
    if param.dgtype =="DGStd"
        refelem = RefElemStd(bnodes, qnodes, fnodes)
        dg = DGStd(PDE.Nstates, refelem, mesh)
    elseif param.dgtype == "DGFluxDiff"
        refelem = RefElemSBP(bnodes, qnodes, fnodes)
        dg = DGFluxDiff(PDE.Nstates, refelem, mesh)
    elseif param.dgtype == "DGArtVisc"
        refelem = RefElemStd(bnodes, qnodes, fnodes)
        dg = DGArtVisc(PDE.Nstates, refelem, mesh)
    elseif param.dgtype == "DGAddRes"
        refelem = RefElemStd(bnodes, qnodes, fnodes)
        dg = DGAddRes(PDE.Nstates, refelem, mesh)
    else
        error("Unknown DG type!")
    end

    # Initialize BCs and physics object
    BChandler = initialize_BCHandler(dg, param)
    physics = PhysProp(PDE, source, BChandler, numflux, tpflux, artviscmodel, rescorrmodel)

    return ic, physics, dg
end

function set_up_and_solve(param::parameters)
    # setup problem
    ic, physics, dg = set_up_problem(param)

    # create error nodes if we need them
    if isnothing(param.enodes) == false
        enodes = make_nodes(param.enodes)
    end

    # Create our time integrator
    method = parse_parameters_timeintegrator(param)

    # Initialize state
    u0 = initialize_states(ic, dg, physics.PDE)

    # Now, we solve
    output = ODE_solver(u0, dg, physics, method) # dictionary containing solution, final time and the time-history of all requested scalar observables
    output["dg"] = dg # we add dg object to solution dict

    # if we want the error, we compute it
    if param.OOAtest
        output["L2error"] = compute_L2error(output["solution"], output["time"], enodes, ic, dg, physics)
        return output
    else
        return output
    end
end

function refinemesh!(param::parameters)
    param.Neldim *= 2
end


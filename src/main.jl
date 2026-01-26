using Plots

function run(param_file)
    param = parse_parameters(param_file)

    if param.OOAtest
        Nrefinements = param.Nrefinements
        L2error = zeros((Nrefinements,)) # store DOFS and Nrefinements
    else
        Nrefinements = 1
    end

    for isolve = 1:Nrefinements
        if param.OOAtest
            sol, L2error[isolve,:] = set_up_and_solve(param)
            refinemesh!(param)
        else
            sol = set_up_and_solve(param)
        end

        N = length(sol)

        if N == 2
            outputnames = ("solution", "time")
        elseif N == 3
            outputnames = ("solution", "time", "entropy")
        end

        for iname in 1:N
            if param.save
            else
                println(name[i] * " :")
                println(sol[i])
            end
        end
    end

    if param.save
    else
        println("L2Error: ")
        println(L2error)
    end

end

function set_up_and_solve(param::parameters)
    # initialize mesh and nodes
    mesh = initialize_mesh(param)
    bnodes = make_nodes(param.bnodes)
    qnodes = make_nodes(param.qnodes)
    fnodes = make_nodes(param.refelem, param.fnodes)
    enodes = make_nodes(param.enodes)

    # initialize ref element and DG object
    if param.pdetype == "LinAdv"
        Nstates = 1
    elseif param.pdetype == "Burgers"
        Nstates = 1
    elseif para.pdetype == "EulerPerfGas"
        Nstates = 2 + param.dim
    end

    if param.dgtype =="DGStd"
        refelem = RefElemStd(bnodes, qnodes, fnodes)
        dg = DGStd(Nstates, refelem, mesh)
    else
        error("Unknown DG type!")
    end

    # Initialize state and BCs
    u0 = initialize_states(dg, param)
    BChandler = initialize_BCHandler(dg, param)

    # Now, we solve
    sol = ODE_solver(u0, BChandler, dg, param) # tuple containing solution, final time and the time-history of all requested scalar observables

    # if we want the error, we compute it
    if param.OOAtest
        error = [dg.DOF, compute_L2error(sol[1], sol[2], enodes, dg, param)]
        return sol, dg.bpts, error
    else
        return sol, dg.bpts
    end
end

function refinemesh!(param::parameters)
    param.Neldim *= 2
end


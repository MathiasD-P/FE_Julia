using FE_Julia
using DelimitedFiles
using Plots

function test_OOA(param, Nrefinements; filename=nothing)
    param.OOAtest=true
    param.calc_entropy=false

    error = zeros((Nrefinements,3))
    tab_title = ["DOF" "L2 Error" "OOA"]
    output = Dict()

    for i in 1:Nrefinements
        try
            output = set_up_and_solve(param)
            error[i,1:2] = [output["dg"].DOF, output["L2error"]]
        catch e
            error[i,1:2] = [NaN, NaN]
            println("careful! Testcase crached.")
            println(e)
        end

        param.Neldim = param.Neldim * 2 #refine mesh

        # We plot the solution for the mesh with max refinement
        if i == Nrefinements
            if param.pdetype == "EulerPerfGas"
                p1 = plot(output["dg"].bpts, output["solution"], xlabel="x", ylabel="u", labels=["rho" "m" "E"])
            else
                p1 = plot(output["dg"].bpts, output["solution"], xlabel="x", ylabel="u")
            end
            display(p1)
        end
    end

    # We plot the convergence of the L2 error vs DOF
    p2 = scatter(log10.(error[:,1]), log10.(error[:,2]), ylabel="log10(L2)", xlabel="log10(DOF)")
    display(p2)

    # Compute OOA and package
    error[2:end,3] .= (log.(error[2:end,2]) .- log.(error[1:end-1,2])) ./ (log.(error[2:end,1]) .- log.(error[1:end-1,1]))

    # Now we print in REPL or save
    if isnothing(filename)
        writedlm(stdout, tab_title)
        writedlm(stdout, error)
        print("\n")
    else
        if !(isdir("outputs"))
            mkdir("outputs")
        end

        open("outputs/" * filename, "w") do IO
            writedlm(IO, tab_title)
            writedlm(IO, error)
            print(IO, "\n")
        end
    end

    return error
end

function test_entropy(param; filename=nothing)
    if !(param.calc_entropy)
        param.calc_entropy=true
    end

    output = set_up_and_solve(param)

    # We plot the solution
    if param.pdetype == "EulerPerfGas"
        p1 = plot(output["dg"].bpts, output["solution"], xlabel="x", ylabel="u", labels=["rho" "m" "E"])
    else
        p1 = plot(output["dg"].bpts, output["solution"], xlabel="x", ylabel="u")
    end
    display(p1)

    # We plot the entropy time series
    p2 = plot(collect(param.dt:param.dt:(param.dt * param.Nsteps)), output["entropy"], xlabel="t", ylabel="S", legend=false)
    display(p2)

    entropy_error = maximum(abs.(output["entropy"] .- output["entropy"][1]))
    entropy = [collect(param.dt:param.dt:(param.dt * param.Nsteps)) output["entropy"]]

    # Now we print in REPL or save
    if isnothing(filename)
        # writedlm(stdout, ["time" "Entropy"])
        # writedlm(stdout, entropy)
        # print("\n")
        print("Entropy Error: ")
        println(entropy_error)

    else
        if !(isdir("outputs"))
            mkdir("outputs")
        end

        open("outputs/" * filename, "w") do IO
            writedlm(IO, ["time" "Entropy"])
            writedlm(IO, entropy)
            print(IO, "\n")
            print(IO, "Entropy Error: ")
            println(IO, entropy_error)
        end
    end

    return entropy_error
end

function test_viscosity(param, Nrefinements; filename=nothing)
    if param.dgtype != "DGArtVisc"
        error("This test can only be applied to DGArtVisc!")
    end
    if isnothing(param.enodes)
        error("You must provide enodes to project ICs.")
    end

    error = zeros((Nrefinements,6))
    tab_title = ["DOF" "Delta" "OOA" "Visc" "OOA" "Den"]

    for i in 1:Nrefinements
        # initialize mesh and nodes
        mesh = FE_Julia.initialize_mesh(param)
        bnodes = FE_Julia.make_nodes(param.bnodes)
        qnodes = FE_Julia.make_nodes(param.qnodes)
        fnodes = FE_Julia.make_nodes(param.refelem, param.fnodes)
        enodes = FE_Julia.make_nodes(param.enodes)

        # initialize ref element and DG object
        if param.pdetype == "LinAdv"
            Nstates = 1
        elseif param.pdetype == "Burgers"
            Nstates = 1
        elseif param.pdetype == "EulerPerfGas"
            Nstates = 2 + param.dim
        end

        refelem = RefElemStd(bnodes, qnodes, fnodes)
        dg = DGArtVisc(Nstates, refelem, mesh)

        # Initialize states on enodes and project
        refelemL2 = RefElemStd(bnodes, enodes, fnodes)
        epts = reduce(vcat, FE_Julia.mapping(dg.mesh, refelemL2.qnodes, ielem) for ielem in 1:dg.mesh.Nel)

        # Initialize state and BCs
        u0 = FE_Julia.block_matmul(refelemL2.Ph, FE_Julia.initialize_states(dg, param, epts), mesh.Nel) # Project initial conditions
        BChandler = FE_Julia.initialize_BCHandler(dg, param)

        # Compute viscosity and entropy deficit
        residual, debug_data = FE_Julia.build_residual!(copy(u0), u0, 0.0, BChandler, dg, param, true)
        error[i,1] = dg.DOF
        error[i,2] = maximum(abs.(debug_data["delta"]))
        error[i,4] = maximum((debug_data["visc"]))
        error[i,6] = minimum((debug_data["den"]))

        param.Neldim = round(Int, param.Neldim * 1.1 + 0.5) #refine mesh
    end

    # We plot the convergence of the L2 error vs DOF
    p2 = scatter(log10.(error[:,1]), log10.(error[:,2]), ylabel="log10(L2)", xlabel="log10(DOF)", label="delta")
    plot!(p2, [log10.(error[1,1]), log10.(error[end,1])], ([log10.(error[1,1]), log10.(error[end,1])] .- log10.(error[end,1])) .* -9 .+ log10.(error[end,2]))
    scatter!(p2, log10.(error[:,1]), log10.(error[:,4]), label="visc")
    plot!(p2, [log10.(error[1,1]), log10.(error[end,1])], ([log10.(error[1,1]), log10.(error[end,1])] .- log10.(error[end,1])) .* -6.5 .+ log10.(error[end,4]))
    scatter!(p2, log10.(error[:,1]), log10.(error[:,6]), label="den")
    display(p2)

    # Compute OOA and package
    error[2:end,3] .= (log.(error[2:end,2]) .- log.(error[1:end-1,2])) ./ (log.(error[2:end,1]) .- log.(error[1:end-1,1]))
    error[2:end,5] .= (log.(error[2:end,4]) .- log.(error[1:end-1,4])) ./ (log.(error[2:end,1]) .- log.(error[1:end-1,1]))

    # Now we print in REPL or save
    if isnothing(filename)
        writedlm(stdout, tab_title)
        writedlm(stdout, error)
        print("\n")
    else
        if !(isdir("outputs"))
            mkdir("outputs")
        end

        open("outputs/" * filename, "w") do IO
            writedlm(IO, tab_title)
            writedlm(IO, error)
            print(IO, "\n")
        end
    end

    return error[end,3], error[end,5]
end

function test_dtmax(param, Tfinal, tol)
    function you_do_explode(solution::Array)
        check = (any(isnan, solution) || any(x -> abs(x) > param.maxval, solution))
        return check
    end

    Nsolvesmax = 1000

    dtlim = copy(param.dtlim)
    if dtlim[1] >= dtlim[2]
        error("You must initialize dt1 < dt2!")
    end

    # Test if initial dtlims work
    crash_flag = false
    try
        param.dt = dtlim[1]
        param.Nsteps = round(Int, Tfinal / param.dt) + 1
        output = set_up_and_solve(param)
        
        if !you_do_explode(output["solution"])
            println("dt1 passes!")
        else
            crash_flag = true
        end

    catch e
        println(e)
        error("Test case should not crash for dt1")
    end
    if crash_flag
        error("Test case should not crash for dt1")
    end


    crash_flag = true
    try
        param.dt = dtlim[2]
        param.Nsteps = round(Int, Tfinal / param.dt) + 1
        output = set_up_and_solve(param)
        
        if you_do_explode(output["solution"])
            println("dt2 crashes!")
        else
            crash_flag = false
        end

    catch e
        println(e)
        println("dt2 crashes!")
    end
    if !crash_flag
        error("Test case should crash for dt2")
    end

    Nsolves = 0
    while abs(dtlim[2] - dtlim[1]) > tol
        dt = 0.5 * sum(dtlim)
        try
            param.dt = dt
            param.Nsteps = round(Int, Tfinal / param.dt) + 1
            output = set_up_and_solve(param)

            if you_do_explode(output["solution"])
                dtlim[2] = dt
            else
                dtlim[1] = dt
            end
        catch e
            println(e)
            dtlim[2] = dt
        end

        Nsolves +=1
        if Nsolves >= Nsolvesmax
            println("Maximum number of iterations exceeded!")
            break
        end
    end

    return dtlim
end
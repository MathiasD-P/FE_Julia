using FE_Julia
using DelimitedFiles
using Plots

function test_OOA(param, Nrefinements; filename=nothing)
    param.OOAtest=true
    param.calc_entropy=false

    error = zeros((Nrefinements,3))
    tab_title = ["DOF" "L2 Error" "OOA"]

    for i in 1:Nrefinements
        output = set_up_and_solve(param)
        error[i,1:2] = [output["dg"].DOF, output["L2error"]]

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
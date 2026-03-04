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

    return error[end,3]
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
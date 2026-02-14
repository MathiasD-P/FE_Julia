using FE_Julia
using Plots

function test_DGFluxDiff_Burgers_1D()

    param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="(7)-GL",
                        qnodes="(8)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_split",
                        twoptfluxtype="EC_split",
                        ICname="exp_1state",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.00002,
                        calc_entropy=true)


    output = set_up_and_solve(param)

    p = plot(output["dg"].bpts, output["solution"], xlabel="x", ylabel="u", legend=false)
    display(p)

    p = plot(collect(0.00002:0.00002:0.2), output["entropy"], xlabel="t", ylabel="S", legend=false)
    display(p)

    entropy_error = maximum(abs.(output["entropy"] .- output["entropy"][1]))

    println("Entropy Error is:")
    println(entropy_error)

    return entropy_error
end
using FE_Julia
using Plots

function test_DGFluxDiff_Chand_Euler_1D_OOA()

    N = 6

    param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="(3)-GL",
                        qnodes="(4)-GLL",
                        enodes="(10)-GLL",
                        fnodes="(1)-GLL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="EC_Chandrashekar",
                        twoptfluxtype="EC_Chandrashekar",
                        ICname="IsentropicDensityWave",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=1000,
                        dt=0.001,
                        OOAtest=true,
                        gamma=1.4)

    error = zeros((N,2))
    for i in 1:N
        output = set_up_and_solve(param)
        error[i,:] = [output["dg"].DOF, output["L2error"]]
        param.Neldim = param.Neldim * 2

        if i == N
            p2 = plot(output["dg"].bpts, output["solution"], xlabel="x", ylabel="u")
            display(p2)
        end
    end

    p1 = scatter(log10.(error[:,1]), log10.(error[:,2]), ylabel="log10(L2)", xlabel="log10(DOF)")
    display(p1)

    OOA = (log(error[end,2]) - log(error[end-1,2])) / (log(error[end,1]) - log(error[end-1,1]))
    println("Oder of accuracy is:")
    println(OOA)

    return OOA
end
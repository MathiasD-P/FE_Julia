using FE_Julia
using Plots

function test_DGStd_Chand_Euler_1D()

    N = 5

    param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="(5)-GL",
                        qnodes="(5)-GLL",
                        enodes="(10)-GLL",
                        fnodes="(1)-GLL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="EC_Chandrashekar",
                        ICname="IsentropicDensityWave",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        gamma=1.4)

    error = zeros((N,2))
    for i in 1:N
        sol, pts, error[i,:] = set_up_and_solve(param)
        param.Neldim = param.Neldim * 2

        if i == N
            p2 = plot(pts, sol[1], xlabel="x", ylabel="u")
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
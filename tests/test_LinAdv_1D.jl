using FE_Julia
using Plots

function test_LinAdv_1D()

    N = 6

    param = parameters(
                        pdetype="LinAdv",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="(5)-GL",
                        qnodes="(6)-GL",
                        enodes="(10)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="upwind",
                        ICname="exp_1state",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        a=1.2)

    error = zeros((N,2))
    for i in 1:N
        sol, pts, error[i,:] = set_up_and_solve(param)
        param.Neldim = param.Neldim * 2

        if i == N
            p2 = plot(pts, sol[1])
            display(p2)
        end
    end

    p1 = scatter(log10.(error[:,1]), log10.(error[:,2]))
    display(p1)

    println("Oder of accuracy is:")
    println((log(error[end,2]) - log(error[end-1,2])) / (log(error[end,1]) - log(error[end-1,1])))

    return
end
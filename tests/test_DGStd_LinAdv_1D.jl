using FE_Julia
using Plots

function test_LinAdv_1D()

    N = 6

    param = parameters(
                        pdetype="LinAdv",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="(6)-GL",
                        qnodes="(7)-GL",
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
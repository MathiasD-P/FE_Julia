using FE_Julia

using Plots
using LegendrePolynomials
using LinearAlgebra
using LaTeXStrings

# Default plotting parameters
default(
    fontfamily = "Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(10, "Computer Modern"),
    tickfont = font(8, "Computer Modern"),
    legendfont = font(8, "Computer Modern")
)

function make_parameters(myDG, p, N, node_name)

    param = parameters(
                            pdetype="Burgers",
                            dim=1,
                            bnodes="(" * string(p+1) * ")-"* node_name,
                            qnodes="(" * string(p+1) * ")-"* node_name,
                            enodes="(12)-GL",
                            fnodes="(1)-GL",
                            refelem="interval",
                            domain="unit_interval_linear",
                            Neldim=N,
                            numfluxtype="LF",
                            ICname="sin_1state",
                            k=2*pi,
                            phi=0.1,
                            av=0.01,
                            BCname="periodic",
                            ODE_solver="LSERK45",
                            Nsteps = 1000,
                            OOAtest=true)
    
    param.dt = 0.6 * (1.0 / (2*pi)) / param.Nsteps

    if myDG == "DGExact"
        param.dgtype = "DGStd"
        param.qnodes = "(" * string(2*p+1) * ")-GL" # overwrite nodes
    else
        param.dgtype = myDG

        if myDG == "DGFluxDiff"
            param.twoptfluxtype = "EC_split"
        elseif myDG == "DGArtVisc"
            param.AVcoeff = "AVdissip"
        elseif myDG == "DGAddRes"
            param.Rescorr = "RescorrEC"
        end
    end

    return param
end

function research_test(p, Nrefinements, node_name, DGnames)
    Nel = 2 # initialization for refinement
    factor = 2 # factor for mesh refinement
    interrors = Dict(myDG => zeros((Nrefinements, 3)) for myDG in DGnames) # initialize dictionary
    interrors["DGExact"] = zeros((Nrefinements, 3))

    for iref = 1:Nrefinements

        paramExact = make_parameters("DGExact", p, Nel, node_name)

        params = Dict()
        for myDG in DGnames
            params[myDG] = make_parameters(myDG, p, Nel, node_name)
        end

        for myDG in DGnames
            output = set_up_and_solve(params[myDG])

            interrors[myDG][iref,2] = output["L2error"]
            interrors[myDG][iref,1] = output["dg"].DOF
            
            if iref != 1
                interrors[myDG][iref, 3] = log(interrors[myDG][iref,2] / interrors[myDG][iref-1,2]) / log(interrors[myDG][iref,1] / interrors[myDG][iref-1,1])
            end

            println(myDG)
        end

        output = set_up_and_solve(paramExact)
        interrors["DGExact"][iref,2] = output["L2error"]
        interrors["DGExact"][iref,1] = output["dg"].DOF
        if iref != 1
            interrors["DGExact"][iref, 3] = log(interrors["DGExact"][iref,2] / interrors["DGExact"][iref-1,2]) / log(interrors["DGExact"][iref,1] / interrors["DGExact"][iref-1,1])
        end

        Nel *= factor
    end

    return interrors
end

function plot_results(interrors, DGnames, p)
    plt_err = plot(
                xscale=:log10,
                yscale=:log10,
                xlabel="DOF",
               ylabel="L2 Error",
               legend=:topright,
               dpi = 700)
    
    for myDG in DGnames

        if myDG == "DGStd"
            mylabel = "Standard DG"
            mymarker = :circle
            mycolor = :black
            mymarkersize = 5
        elseif myDG == "DGFluxDiff"
            mylabel = "Flux Differencing"
            mymarker = :utriangle
            mycolor = :red
            mymarkersize = 5
        elseif myDG == "DGArtVisc"
            mylabel = "Artificial Viscosity"
            mymarker = :star5
            mycolor = :blue
            mymarkersize = 5
        elseif myDG == "DGAddRes"
            mylabel = "Dealiased Split"
            mymarker = :diamond
            mycolor = :green
            mymarkersize = 5
        end
        scatter!(plt_err, interrors[myDG][:,1], interrors[myDG][:,2], label=mylabel, markershape=mymarker, markersize=mymarkersize, color=mycolor)
    end

    scatter!(plt_err, interrors["DGExact"][:,1], interrors["DGExact"][:,2], label="Exact Integration", markershape=:rect, markersize=3.5, color=:orange)

    plot!(plt_err, interrors["DGStd"][end-2:end,1], 2 .* interrors["DGStd"][end-2,2] .* (interrors["DGStd"][end-2:end,1] ./ interrors["DGStd"][end-2,1]).^(-(p+1)), label=L"p+1", color=:black, linestyle=:dash)

    return plt_err
end

DGnames = ["DGStd", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
p = 5
Nrefinements = 8
node_name = "GLL"

errors = research_test(p, Nrefinements, node_name, DGnames)
plt = plot_results(errors, DGnames, p)
# savefig(plt, "L2error" * node_name * "p" * string(p) * ".pdf")
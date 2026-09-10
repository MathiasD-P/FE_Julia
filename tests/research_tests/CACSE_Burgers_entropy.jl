using FE_Julia
using Plots
using LaTeXStrings

include("../validation_tests_parameters.jl")
include("../validation_tests.jl")

# Default plotting parameters
default(
    fontfamily="Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(12, "Computer Modern"),
    tickfont = font(10, "Computer Modern"),
    legendfont = font(10, "Computer Modern")
)

function make_parameters(bnodes, qnodes)
    param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        bnodes=bnodes,
                        qnodes=qnodes,
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_split",
                        twoptfluxtype="EC_split",
                        AVcoeff="AVEC",
                        Rescorr="RescorrEC",
                        ICname="exp_1state",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.00002,
                        calc_entropy=true)

    return param
end

function research_test(bnodes, qnodes, DGnames, colors, filename=nothing)
    # Construct parameters
    param = make_parameters(bnodes, qnodes)

    # Testing
    results = Dict()
    for myDG in DGnames
        param.dgtype = myDG
        results[myDG] = set_up_and_solve(param)
    end

    # plotting

    plt_sol = plot(
               xlabel=L"x",
               ylabel=L"u",
               dpi=500)

    for (colori, myDG) in enumerate(DGnames)
        plot!(plt_sol, results[myDG]["dg"].bpts, results[myDG]["solution"], color=colors[colori])
    end
    display(plt_sol)

    plt_ent = plot(
               xlabel=L"t",
               ylabel=L"U(t)-U(0)",
               title= "PDE: " * param.pdetype * "; Basis: " * bnodes * "; Quad: " * qnodes,
               dpi=500)

    for (colori, myDG) in enumerate(DGnames)
        if myDG == "DGArtVisc"
            mylabel = "Artificial Viscosity"
        elseif myDG == "DGAddRes"
            mylabel = "Residual Correction"
        elseif myDG == "DGStd"
            mylabel = "DG"
        elseif myDG == "DGFluxDiff"
            mylabel = "Flux Differencing"
        end
        plot!(plt_ent, collect(param.dt:param.dt:(param.Nsteps*param.dt)), results[myDG]["entropy"] .- results[myDG]["entropy"][1], color=colors[colori], label=mylabel)
    end

    if isnothing(filename)
        display(plt_ent)
    else
        savefig(plt_ent, filename)
    end

    return results
end

bnodes = "(5)-GLL"
qnodes = "(5)-GLL"
DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
colors = [:red, :blue, :green]

resultsGLL= research_test(bnodes, qnodes, DGnames, colors, "outputs/GLL_ent.png")

bnodes = "(5)-GL"
qnodes = "(5)-GL"
DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
colors = [:red, :blue, :green]

resultsGL= research_test(bnodes, qnodes, DGnames, colors, "outputs/GL_ent.png")
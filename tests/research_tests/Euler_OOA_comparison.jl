#####################################################################
# Basically just a plotting code for the OOA for Burgers
#####################################################################

using FE_Julia
using Plots
using LaTeXStrings

include("../validation_tests_parameters.jl")
include("../validation_tests.jl")

# Default plotting parameters
default(
    fontfamily = "Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(10, "Computer Modern"),
    tickfont = font(8, "Computer Modern"),
    legendfont = font(8, "Computer Modern")
)

function research_test(bnodes, qnodes, DGnames, colors, Nrefinements)
    # Construct parameters
    params = Dict()
    for myDG in DGnames
        params[myDG] = make_validation_tests_parameters("test_OOA_" * myDG * "_Chand_Euler_1D_" * "b" * bnodes * "_q" * qnodes)
    end

    # Testing starts here
    errors = Dict()
    for myDG in DGnames
        println(myDG)
        errors[myDG] = test_OOA(params[myDG], Nrefinements)
    end

    # Now, we plot
    plt = plot(xscale = :log10,
               yscale = :log10,
               title= "PDE: " * collect(values(params))[1].pdetype * "; Basis: " * bnodes * "; Quad: " * qnodes,
               xlabel="DOF",
               ylabel="L2 Error")

    for (colori, myDG) in enumerate(DGnames)
        if myDG == "DGArtVisc"
            mylabel = myDG * " (" * params[myDG].AVcoeff * ")"
        elseif myDG == "DGAddRes"
            mylabel = myDG * " (" * params[myDG].Rescorr * ")"
        else
            mylabel = myDG
        end
        scatter!(plt, (errors[myDG][:,1]), (errors[myDG][:,2]), color=colors[colori], label=mylabel)
    end

    display(plt)

    # We compare error
    myplots = []
    for iDOF in 1:Nrefinements
        plt = plot(legend=false, xticks=false)
        for (colori, myDG) in enumerate(DGnames)
            scatter!(plt, [errors[myDG][iDOF,1]], [errors[myDG][iDOF,2]], color=colors[colori], xlims = (errors[myDG][iDOF,1]-0.01, errors[myDG][iDOF,1]+0.01))
        end
        push!(myplots, plt)
    end
    plt=plot(myplots..., layout=(1,length(myplots)))
    display(plt)
end

# Testing
bnodes = "(4)-GLL"
qnodes = "(4)-GLL"
Nrefinements = 7
DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes", "DGStd"]
colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442"]

errors = research_test(bnodes, qnodes, DGnames, colors, Nrefinements)

# # Testing
# bnodes = "(4)-GLL"
# qnodes = "(4)-GL"
# Nrefinements = 7
# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes", "DGStd"]
# colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442"]

# errors = research_test(bnodes, qnodes, DGnames, colors, Nrefinements)
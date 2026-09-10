#####################################################################
# Basically just a plotting code for the OOA for Burgers
#####################################################################

using FE_Julia
using Plots
using LaTeXStrings
using Printf

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

function research_test(bnodes, qnodes, DGnames, Nrefinements)
    # Construct parameters
    params = Dict()
    for myDG in DGnames
        params[myDG] = make_validation_tests_parameters("test_OOA_" * myDG * "_Burgers_1D_" * "b" * bnodes * "_q" * qnodes)
    end

    # Testing starts here
    errors = Dict()
    for myDG in DGnames
        println(myDG)
        errors[myDG] = test_OOA(params[myDG], Nrefinements)
    end

    return errors
end

function plot_test(errors, bnodes, qnodes, colors, DGnames, filename=nothing)
    # reconstruct parameters
    params = Dict()
    for myDG in DGnames
        params[myDG] = make_validation_tests_parameters("test_OOA_" * myDG * "_Burgers_1D_" * "b" * bnodes * "_q" * qnodes)
    end

    order = parse(Float64, string(bnodes[2]))

    plt = plot(xscale = :log10,
               yscale = :log10,
               title= "PDE: " * collect(values(params))[1].pdetype * "; Basis: " * bnodes * "; Quad: " * qnodes,
               xlabel="DOF",
               ylabel="L2 Error")

    # We compare error (all DOFs)
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
        scatter!(plt, (errors[myDG][:,1]), (errors[myDG][:,2]), color=colors[colori], label=mylabel)
    end
    aDG = collect(keys(errorsGL))[1]
    plot!(plt, (errors[aDG][end-2:end,1]), 2 * errors[aDG][end,2] .* (errors[aDG][end-2:end,1] ./ errors[aDG][end,1]).^-order, color=:black, linestyle=:dash, label="Order " * string(order), ylims=(1e-10,30), dpi=500)
    display(plt)

    if isnothing(filename)
        display(plt)
    else
        savefig(plt, filename * ".png")
    end

    # We compare error (zoom over 3 first DOFs)
    plt = plot(xscale = :log10,
               yscale = :log10,
               title= "PDE: " * collect(values(params))[1].pdetype * "; Basis: " * bnodes * "; Quad: " * qnodes,
               xlabel="DOF",
               ylabel="L2 Error")

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
        scatter!(plt, (errors[myDG][1:3,1]), (errors[myDG][1:3,2]), color=colors[colori], label=mylabel)
    end
    display(plt)
end

function print_table(errors, DGnames)
    DOF = size(errors[DGnames[1]], 1)
    for iDOF in 1:DOF
        s = ""
        for (iDG, myDG) in enumerate(DGnames)
            s *= @sprintf("%.1e", errors[myDG][iDOF,2])

            if iDG != length(DGnames)
                s *= " & "
            else
                s *= " \\\\"
            end
        end

        println(s)
    end
end

# Testing GL
# bnodes = "(5)-GL"
# qnodes = "(5)-GL"
# Nrefinements = 7
# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes", "DGStd"]
# colors = [:red, :blue, :green, :orange]

# errorsGL = research_test(bnodes, qnodes, DGnames, Nrefinements)
# plot_test(errorsGL, bnodes, qnodes, colors, DGnames, "outputs/OOA_GL")

# Testing GLL
bnodes = "(5)-GLL"
qnodes = "(5)-GLL"
Nrefinements = 7
DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes", "DGStd"]
colors = [:red, :blue, :green, :orange]

errorsGLL = research_test(bnodes, qnodes, DGnames, Nrefinements)
plot_test(errorsGLL, bnodes, qnodes, colors, DGnames, "outputs/OOA_GLL")
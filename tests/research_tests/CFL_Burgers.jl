
using FE_Julia
using Plots
using LaTeXStrings

include("../validation_tests.jl")

# Default plotting parameters
default(
    fontfamily = "Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(10, "Computer Modern"),
    tickfont = font(8, "Computer Modern"),
    legendfont = font(8, "Computer Modern")
)

function research_test(param, tol, Tfinal, DGnames, colors, Nrefinements)
    # code execution for each dgname
    out = Dict()
    Neldim_init = param.Neldim
    # dtlim_init = copy(param.dtlim)

    for mydg in DGnames
        out[mydg] = zeros(Nrefinements, 2)

        for i in 1:Nrefinements
            param.dgtype = mydg
            out[mydg][i,:] .= test_dtmax(param, Tfinal, tol)
            param.Neldim *= 2
        end

        # make sure to reset parameters
        param.Neldim = Neldim_init
    end

    # plot results
    dx = 1 / param.Neldim .* 0.5.^collect(range(0,Nrefinements-1, Nrefinements))

    plt = plot(title= "PDE: " * param.pdetype * "; Basis: " * param.bnodes * "; Quad: " * param.qnodes,
               xlabel="dx",
               ylabel="dt/dx",
               xscale = :log10,
               legend = :bottomleft)

    for (colori, myDG) in enumerate(DGnames)
        if myDG == "DGArtVisc"
            mylabel = myDG * " (" * param.AVcoeff * ")"
        elseif myDG == "DGAddRes"
            mylabel = myDG * " (" * param.Rescorr * ")"
        else
            mylabel = myDG
        end

        scatter!(plt, dx, 0.5 .* (out[myDG][:,1] .+ out[myDG][:,2]) ./ dx, yerr=abs.(0.5 .* (out[myDG][:,1] .- out[myDG][:,2])) ./ dx, color=colors[colori], label=mylabel)
    end

    display(plt)

    return dx, out
end

# (5)-GLL quad ; Gaussian

# param = parameters(                  
#                     pdetype="Burgers",
#                     dgtype=nothing,
#                     dim=1,
#                     bnodes="(5)-GLL",
#                     qnodes="(5)-GLL",
#                     fnodes="(1)-GL",
#                     refelem="interval",
#                     domain="unit_interval_linear",
#                     Neldim=10,
#                     numfluxtype="LF",
#                     twoptfluxtype="EC_split",
#                     Rescorr="RescorrEC",
#                     AVcoeff="AVdissip",
#                     ICname="exp_1state",
#                     BCname="periodic",
#                     ODE_solver="LSERK45",
#                     maxval=50,
#                     dtlim=[0.0005, 0.1])

# tol = 0.0001
# Tfinal = 5.0
# Nrefinements = 4
# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes", "DGStd"]
# colors = [:red, :blue, :green, :orange]

# dx, dt = research_test(param, tol, Tfinal, DGnames, colors, Nrefinements)

# # (5)-GL quad ; Gaussian

# param = parameters(                  
#                     pdetype="Burgers",
#                     dgtype=nothing,
#                     dim=1,
#                     bnodes="(5)-GLL",
#                     qnodes="(5)-GL",
#                     fnodes="(1)-GL",
#                     refelem="interval",
#                     domain="unit_interval_linear",
#                     Neldim=10,
#                     numfluxtype="LF",
#                     twoptfluxtype="EC_split",
#                     Rescorr="RescorrEC",
#                     AVcoeff="AVdissip",
#                     ICname="exp_1state",
#                     BCname="periodic",
#                     ODE_solver="LSERK45",
#                     maxval=50,
#                     dtlim=[0.0005, 0.1])

# tol = 0.0001
# Tfinal = 5.0
# Nrefinements = 4
# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes", "DGStd"]
# colors = [:red, :blue, :green, :orange]

# dx, dt = research_test(param, tol, Tfinal, DGnames, colors, Nrefinements)

# (5)-GLL quad ; sin_1state

param = parameters(                  
                    pdetype="Burgers",
                    dgtype=nothing,
                    dim=1,
                    bnodes="(5)-GLL",
                    qnodes="(5)-GLL",
                    fnodes="(1)-GL",
                    refelem="interval",
                    domain="unit_interval_linear",
                    Neldim=5,
                    numfluxtype="LF",
                    twoptfluxtype="EC_split",
                    Rescorr="RescorrEC",
                    AVcoeff="AVdissip",
                    ICname="sin_1state",
                    BCname="periodic",
                    ODE_solver="LSERK45",
                    maxval=50,
                    dtlim=[0.0001, 0.1],
                    k=[4.0])

tol = 0.0001
Tfinal = 5.0
Nrefinements = 4
DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
colors = [:red, :blue, :green, :orange]

dx, dt = research_test(param, tol, Tfinal, DGnames, colors, Nrefinements)
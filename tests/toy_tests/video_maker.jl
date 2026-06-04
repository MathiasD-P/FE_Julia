using FE_Julia
using LaTeXStrings
using Plots

default(
    fontfamily="Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(12, "Computer Modern"),
    tickfont = font(10, "Computer Modern"),
    legendfont = font(10, "Computer Modern")
)

function make_video(param::parameters, Nsave::Integer, chiplt::Matrix, pltrefpts::Array, mystyle;  filename=nothing) # CAREFUL, VIDEO MAKER ONLY WORKS FOR AUTONOMOUS SYSTEMS AT THE MOMENT
    # Change default Nsteps to compute until next frame
    Nsteps_bysave = -floor(Int, -param.Nsteps / Nsave)
    param.Nsteps = Nsteps_bysave

    # Initialize problem
    u0, BChandler, dg = set_up_problem(param)
    t = 0

    # Compute plotting points
    pltpts = reduce(vcat, FE_Julia.mapping(dg.mesh, pltrefpts, ielem) for ielem in 1:dg.mesh.Nel)

    for isave = 1:(Nsave+1)
        if isave != 1
            output = ODE_solver(u0, BChandler, dg, param)
            u0 .= output["solution"]
            t += output["time"]
        end

        uplt = FE_Julia.block_matmul(chiplt, u0, dg.mesh.Nel)
        plt = plot_snapshot(pltpts, uplt, t, param, mystyle)
        
        if isnothing(filename)
            display(plt)
        else
            savefig(plt, filename * "/u" * string(isave) * ".png")
        end
    end

    return
end

function plot_snapshot(pltpts, u0, t, param, mystyle)
    if (param.domain == "unit_interval_linear") & (size(u0, 2) == 1) # 1D and 1 state
        plt = plot(pltpts, u0 ; mystyle..., title = LaTeXString("t = " * string(round(t, digits=4))), dpi=500)
    end

    return plt
end

# Videos for C3SE2026
# # Burgulence

param = parameters(
                    pdetype="Burgers",
                    dgtype = "DGStd",
                    dim=1,
                    bnodes="(5)-GL",
                    qnodes="(5)-GL",
                    fnodes="(1)-GL",
                    refelem="interval",
                    domain="unit_interval_linear",
                    Neldim=50,
                    numfluxtype="LF",
                    ICname="Burgulence",
                    BCname="periodic",
                    ODE_solver="LSERK45",
                    Nsteps=3000,
                    dt=0.0001,
                    kmax=50)

pltrefpts = FE_Julia.evaluate(make_nodes("(25)-GLL"))
chiplt = FE_Julia.Lagrange_Vandermonde1D(pltrefpts, FE_Julia.evaluate(make_nodes(param.bnodes)))
mystyle = (xlabel=L"x", ylabel=L"u", color=:red, legend=nothing, xlims=(-0.5,0.5), ylims=(-0.1,2.5), size=(700, 400))

if !(isdir("outputs/video_Burgulence"))
    mkdir("outputs/video_Burgulence")
end
make_video(param, 30, chiplt, pltrefpts, mystyle, filename="outputs/video_Burgulence")


# Levesque pathology

param = parameters(
                    pdetype="Burgers",
                    dgtype = "DGStd",
                    dim=1,
                    bnodes="(5)-GL",
                    qnodes="(5)-GL",
                    fnodes="(1)-GL",
                    refelem="interval",
                    domain="unit_interval_linear",
                    Neldim=50,
                    numfluxtype="upwind",
                    ICname="rarefaction_1state",
                    BCname="unit_rarefaction",
                    ODE_solver="LSERK45",
                    Nsteps=2000,
                    dt=0.0001)

pltrefpts = FE_Julia.evaluate(make_nodes("(25)-GLL"))
chiplt = FE_Julia.Lagrange_Vandermonde1D(pltrefpts, FE_Julia.evaluate(make_nodes(param.bnodes)))
mystyle = (xlabel=L"x", ylabel=L"u", color=:red, label="Standard DG with upwind flux", xlims=(-0.5,0.5), ylims=(-1.1,1.1), size=(700, 400))

if !(isdir("outputs/video_Levesque_path"))
    mkdir("outputs/video_Levesque_path")
end
make_video(param, 20, chiplt, pltrefpts, mystyle, filename="outputs/video_Levesque_path")

param = parameters(
                    pdetype="Burgers",
                    dgtype = "DGFluxDiff",
                    dim=1,
                    bnodes="(5)-GL",
                    qnodes="(5)-GL",
                    fnodes="(1)-GL",
                    refelem="interval",
                    domain="unit_interval_linear",
                    Neldim=50,
                    numfluxtype="LF",
                    twoptfluxtype="EC_split",
                    Rescorr="RescorrEC",
                    ICname="rarefaction_1state",
                    BCname="unit_rarefaction",
                    ODE_solver="LSERK45",
                    Nsteps=2000,
                    dt=0.0001)

pltrefpts = FE_Julia.evaluate(make_nodes("(25)-GLL"))
chiplt = FE_Julia.Lagrange_Vandermonde1D(pltrefpts, FE_Julia.evaluate(make_nodes(param.bnodes)))
mystyle = (xlabel=L"x", ylabel=L"u", color="#009E73", label="Entropy Stable DG", xlims=(-0.5,0.5), ylims=(-1.1,1.1), size=(700, 400))

if !(isdir("outputs/video_Levesque_ES"))
    mkdir("outputs/video_Levesque_ES")
end
make_video(param, 20, chiplt, pltrefpts, mystyle, filename="outputs/video_Levesque_ES")

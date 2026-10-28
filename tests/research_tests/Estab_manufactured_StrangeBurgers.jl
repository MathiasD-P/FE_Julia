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

function make_parameters(myDG, p, node_name, Nel)

    param = parameters(
                            pdetype="Burgers",
                            dim=1,
                            bnodes="(" * string(p+1) * ")-"* node_name,
                            qnodes="(" * string(p+1) * ")-"* node_name,
                            qsnodes="(" * string(3*p+1) * ")-GL",
                            enodes="(" * string(3*p+1) * ")-GL",
                            fnodes="(1)-GL",
                            refelem="interval",
                            domain="unit_interval_linear",
                            Neldim=Nel,
                            numfluxtype="LF",
                            BCname="periodic",
                            ODE_solver="LSERK45",
                            ICname="StrangeBurgers",
                            sourcename="StrangeBurgers",
                            Nsteps = 30000,
                            OOAtest=true)
    
    param.dt = 3.0 / param.Nsteps

    if node_name == "GLL"
        param.strangecoeff = [1.1, 2.0, 2*pi, 4*pi, 1.0, 1.0, 0.01]
    elseif node_name == "GL"
        param.strangecoeff = [1.1, 2.9, 4*pi, 6*pi, 1.0, 1.0, 0.01]
        param.strangecoeff = [1.1, 2.0, 2*pi, 4*pi, 1.0, 1.0, 0.01]
    end

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

function compute_error(DGnames, p, node_name, Nrange)
    params = Dict()
    outputs = Dict()
    for myDG in DGnames
        params[myDG] = make_parameters(myDG, p, node_name, 2)
        outputs[myDG] = zeros(length(Nrange), 3)
        outputs[myDG][:,1] .= Nrange .* (p+1)
    end

    for myDG in DGnames
        for (i, Nel) in enumerate(Nrange)
            params[myDG].Neldim = Nel
            myoutput = set_up_and_solve(params[myDG])
            outputs[myDG][i, 2] = myoutput["L2error"]

            if i != 1
                outputs[myDG][i, 3] = (log(outputs[myDG][i, 2]) - log(outputs[myDG][i-1, 2])) / (log(Nrange[i]) - log(Nrange[i-1]))
            end
        end
    end
    return outputs
end

function plot_error(DGnames, outputs)
    plt = plot(xlabel="DOF",
               ylabel="L2 Error",
               xaxis=:log10,
               yaxis=:log10,
               dpi = 500)

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
            mylabel = "Residual Correction"
            mymarker = :diamond
            mycolor = :green
            mymarkersize = 5
        elseif myDG == "DGExact"
            mylabel = "Exact Integration"
            mymarker = :square
            mycolor = :orange
            mymarkersize = 3
        end

        scatter!(plt, outputs[myDG][:,1], outputs[myDG][:,2], label=mylabel, color=mycolor, marker=mymarker, markersize=mymarkersize)
        plot!(plt, outputs[myDG][:,1], outputs[myDG][:,2], linestyle=:dash, color=mycolor, label=false)
    end

    return plt
end

function plot_solutions(DGnames, p, node_name, Nel)
    Delta = 0.15

    params = Dict()
    for myDG in DGnames
        params[myDG] = make_parameters(myDG, p, node_name, Nel)
    end

    # Testing starts here
    outputs = Dict()
    for myDG in DGnames
        println(myDG)
        outputs[myDG] = set_up_and_solve(params[myDG])
    end

    plt_sol = plot(xlabel=L"x",
               ylabel=L"u",
               legend=:topright,
               dpi = 500)

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
            mylabel = "Residual Correction"
            mymarker = :diamond
            mycolor = :green
            mymarkersize = 5
        elseif myDG == "DGExact"
            mylabel = "Exact Integration"
            mymarker = :square
            mycolor = :orange
            mymarkersize = 3
        end

        # We use error nodes for plotting
        chie, _, epts = FE_Julia.extract_volume_quadrature(outputs[myDG]["dg"].refelem.bnodestype, make_nodes(params[myDG].enodes))
        epts = reduce(vcat, FE_Julia.mapping(outputs[myDG]["dg"].mesh, epts, ielem) for ielem in 1:outputs[myDG]["dg"].mesh.Nel)
        solution = FE_Julia.block_matmul(chie, outputs[myDG]["solution"], outputs[myDG]["dg"].mesh.Nel)

        ds = sqrt.((epts[2:end] .- epts[1:end-1]).^2 .+ (solution[2:end] .- solution[1:end-1]).^2)
        ds = [0; ds]
        s = [sum(ds[1:i]) for i in 1:length(ds)]
        s[abs.(s) .> 1e10] .= 0
        markerslicer = unique(i -> round.(Int, s ./ Delta)[i], eachindex(s))

        plot!(plt_sol, epts, solution, label=false, color=mycolor)
        scatter!(plt_sol, epts[markerslicer], solution[markerslicer], label=mylabel, color=mycolor, markershape=mymarker, markersize=mymarkersize, markerstrokewidth = 0)
    end

    # We also plot the exact solution
    x = collect(-0.5:0.01:0.5)
    plot!(plt_sol, x, params["DGExact"].strangecoeff[1] .* sin.(params["DGExact"].strangecoeff[2] .* sin.(params["DGExact"].strangecoeff[3] .* x) .+ params["DGExact"].strangecoeff[4] .* x .+ params["DGExact"].strangecoeff[5] .- 3.0 .* params["DGExact"].strangecoeff[6]) .+ params["DGExact"].strangecoeff[7], color=:black, linestyle=:dot, label="Exact Solution")

    return plt_sol
end

# DGnames = ["DGExact", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=5
# node_name = "GLL"
# Nrange = [12, 14, 16, 18, 20, 40, 80, 160]

# outputsGLL = compute_error(DGnames, p, node_name, Nrange)
# pltGLL = plot_error(DGnames, outputsGLL)
# display(pltGLL)
# savefig(pltGLL, "outputs/L2error_StrangeBurgers" * node_name * "p" * string(p) * ".pdf")

# DGnames = ["DGExact", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=5
# node_name = "GL"
# Nrange = [12, 14, 16, 18, 20, 40, 80, 160]

# outputsGL = compute_error(DGnames, p, node_name, Nrange)
# pltGL = plot_error(DGnames, outputsGL)
# display(pltGL)
# savefig(pltGL, "outputs/L2error_StrangeBurgers" * node_name * "p" * string(p) * ".pdf")

DGnames = ["DGExact", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
p=4
node_name = "GLL"
Nrange = [12, 14, 16, 18, 20, 40, 80, 160]

outputsGLL = compute_error(DGnames, p, node_name, Nrange)
pltGLL = plot_error(DGnames, outputsGLL)
display(pltGLL)
savefig(pltGLL, "outputs/L2error_StrangeBurgers" * node_name * "p" * string(p) * ".pdf")

DGnames = ["DGExact", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
p=4
node_name = "GL"
Nrange = [12, 14, 16, 18, 20, 40, 80, 160]

outputsGL = compute_error(DGnames, p, node_name, Nrange)
pltGL = plot_error(DGnames, outputsGL)
display(pltGL)
savefig(pltGL, "outputs/L2error_StrangeBurgers" * node_name * "p" * string(p) * ".pdf")

# DGnames = ["DGExact", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=5
# node_name = "GLL"
# Nel = 16
# plt = plot_solutions(DGnames, p, node_name, Nel)
# display(plt)
# savefig(pltGL, "outputs/solution_StrangeBurgers" * node_name * "p" * string(p) * ".pdf")

# DGnames = ["DGExact", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=5
# node_name = "GL"
# Nel = 16
# plt = plot_solutions(DGnames, p, node_name, Nel)
# display(plt)
# savefig(pltGL, "outputs/solution_StrangeBurgers" * node_name * "p" * string(p) * ".pdf")

# DGnames = ["DGStd", "DGExact"]
# p=8
# node_name = "GLL"
# Nel = 40
# plt = plot_solutions(DGnames, p, node_name, Nel)

# display(plt)

# For GLL"

# Integrate until t = 3.0

# param = make_parameters("DGStd", 5, "GLL", 16)

# ic, physics, dg = set_up_problem(param)
# physics.source = FE_Julia.StrangeBurgersSource(1.1, 2.0, nothing, 4*pi, nothing, nothing, 0.01)
# ic = FE_Julia.StrangeBurgers(physics.source)
# method = FE_Julia.parse_parameters_timeintegrator(param)
# u0 = initialize_states(ic, dg, physics.PDE)
# u = copy(u0)

# output = ODE_solver(u, dg, physics, method)

# plt = plot(dg.bpts, output["solution"])
# plot!(plt, dg.bpts, u0)
# display(plt)

# For GL

# Integrate until t = 3.0

# param = make_parameters("DGArtVisc", 8, "GL", 10)

# ic, physics, dg = set_up_problem(param)
# physics.source = FE_Julia.StrangeBurgersSource(1.1, 2.8, 4*pi, 6*pi, nothing, 1.0, 0.01)
# ic = FE_Julia.StrangeBurgers(physics.source)
# method = FE_Julia.parse_parameters_timeintegrator(param)
# u0 = initialize_states(ic, dg, physics.PDE)
# u = copy(u0)

# output = ODE_solver(u, dg, physics, method)

# plt = plot(dg.bpts, output["solution"])
# plot!(plt, dg.bpts, u0)
# display(plt)

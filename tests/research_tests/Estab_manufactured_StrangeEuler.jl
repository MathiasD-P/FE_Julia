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
                            pdetype="EulerPerfGas",
                            dim=1,
                            gamma=1.4,
                            bnodes="(" * string(p+1) * ")-"* node_name,
                            qnodes="(" * string(p+1) * ")-"* node_name,
                            qsnodes="(" * string(3*p+1) * ")-GL",
                            enodes="(" * string(3*p+1) * ")-GL",
                            fnodes="(1)-GL",
                            refelem="interval",
                            domain="unit_interval_linear",
                            Neldim=Nel,
                            numfluxtype="ES_Chandrashekar_dissip",
                            ICname = "StrangeEuler",
                            sourcename = "StrangeEuler",
                            BCname="periodic",
                            ODE_solver="LSERK45",
                            Nsteps = 4000,
                            OOAtest = true)
    
    param.dt = 1.0 / param.Nsteps
    
    if node_name == "GLL"
        param.strangecoeff = [1.1, 2.0, 2*pi, 4*pi, 1.0, 1.0, 0.01]
    elseif node_name == "GL"
        param.strangecoeff = [1.1, 2.0, 2*pi, 4*pi, 1.0, 1.0, 0.01]
    end

    if myDG == "DGExact"
        param.dgtype = "DGStdP"
        param.qpnodes = param.qnodes
        param.qnodes = "(" * string(10*p+1) * ")-GL" # overwrite nodes
    else
        param.dgtype = myDG

        if myDG == "DGStdP"
            param.qpnodes = param.qnodes
        elseif myDG == "DGFluxDiff"
            param.twoptfluxtype = "EC_Chandrashekar"
        elseif myDG == "DGArtVisc"
            param.AVcoeff = "AVdissip"
        elseif myDG == "DGAddRes"
            param.Rescorr = "RescorrEC"
        end
    end

    return param
end

function exact_solution(x, t)
    solution = Matrix{Float64}(undef, length(x), 3)

    solution[:,1] .= 1.0
    solution[:,2] .= 1.1 .* sin.(2.0.*sin.(2*pi .* x) .+ 4*pi .* x .+ 1.0 .- 1.0 * t) .+ 0.01
    solution[:,3] .= 1.0 / (1.4 - 1) .+ 0.5 .* solution[:,2].^2

    return solution
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

function plot_error(DGnames, p, outputs)
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

    plot!(plt, outputs["DGAddRes"][end-1:end,1], 2 .* outputs["DGAddRes"][end-1,2] .* (outputs["DGAddRes"][end-1:end,1] ./ outputs["DGAddRes"][end-1,1]).^(-(p+1)), label=L"p+1", color=:black, linestyle=:dash)

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

    plt_list = []

    for isol = 1:3
        if isol == 1
            myylabel = L"\rho"
        elseif isol == 2
            myylabel = L"\rho u"
        elseif isol == 3
            myylabel = L"E"
        end

        plt_sol = plot(xlabel=L"x",
                ylabel=myylabel,
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
            solution = FE_Julia.block_matmul(chie, outputs[myDG]["solution"][:,isol], outputs[myDG]["dg"].mesh.Nel) # we only select momentum

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
        plot!(plt_sol, x, exact_solution(x, 1.0)[:, isol], color=:black, linestyle=:dot, label="Exact Solution")

        push!(plt_list, plt_sol)
    end

    return plt_list
end

# Error plots

# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=5
# node_name = "GLL"
# Nrange = [22, 23, 24, 25, 26, 27, 40, 80, 160]

# outputs = compute_error(DGnames, p, node_name, Nrange)

# pltGLL=plot_error(DGnames, p, outputs)
# display(pltGLL)
# savefig(pltGLL, "outputs/L2error_StrangeEuler" * node_name * "p" * string(p) * ".pdf")

# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=5
# node_name = "GL"
# Nrange = [22, 23, 24, 25, 26, 27, 40, 80, 160]

# outputs = compute_error(DGnames, p, node_name, Nrange)

# pltGL=plot_error(DGnames, p, outputs)
# display(pltGL)
# savefig(pltGL, "outputs/L2error_StrangeEuler" * node_name * "p" * string(p) * ".pdf")

# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=4
# node_name = "GLL"
# Nrange = [22, 23, 24, 25, 26, 27, 40, 80, 160]

# outputs = compute_error(DGnames, p, node_name, Nrange)

# pltGLL=plot_error(DGnames, p, outputs)
# display(pltGLL)
# savefig(pltGLL, "outputs/L2error_StrangeEuler" * node_name * "p" * string(p) * ".pdf")

# DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
# p=4
# node_name = "GL"
# Nrange = [22, 23, 24, 25, 26, 27, 40, 80, 160]

# outputs = compute_error(DGnames, p, node_name, Nrange)

# pltGL=plot_error(DGnames, p, outputs)
# display(pltGL)
# savefig(pltGL, "outputs/L2error_StrangeEuler" * node_name * "p" * string(p) * ".pdf")

# Solutions

DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
p=4
node_name = "GLL"
Nel = 22
pltlist = plot_solutions(DGnames, p, node_name, Nel)
display(pltlist[1])
display(pltlist[2])
display(pltlist[3])

savefig(pltlist[1], "outputs/solution_StrangeEuler_density" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[2], "outputs/solution_StrangeEuler_momentum" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[3], "outputs/solution_StrangeEuler_energy" * node_name * "p" * string(p) * ".pdf")

DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
p=4
node_name = "GL"
Nel = 22
pltlist = plot_solutions(DGnames, p, node_name, Nel)
display(pltlist[1])
display(pltlist[2])
display(pltlist[3])

savefig(pltlist[1], "outputs/solution_StrangeEuler_density" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[2], "outputs/solution_StrangeEuler_momentum" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[3], "outputs/solution_StrangeEuler_energy" * node_name * "p" * string(p) * ".pdf")

DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
p=5
node_name = "GLL"
Nel = 16
pltlist = plot_solutions(DGnames, p, node_name, Nel)
display(pltlist[1])
display(pltlist[2])
display(pltlist[3])

savefig(pltlist[1], "outputs/solution_StrangeEuler_density" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[2], "outputs/solution_StrangeEuler_momentum" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[3], "outputs/solution_StrangeEuler_energy" * node_name * "p" * string(p) * ".pdf")

DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
p=5
node_name = "GL"
Nel = 16
pltlist = plot_solutions(DGnames, p, node_name, Nel)
display(pltlist[1])
display(pltlist[2])
display(pltlist[3])

savefig(pltlist[1], "outputs/solution_StrangeEuler_density" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[2], "outputs/solution_StrangeEuler_momentum" * node_name * "p" * string(p) * ".pdf")
savefig(pltlist[3], "outputs/solution_StrangeEuler_energy" * node_name * "p" * string(p) * ".pdf")
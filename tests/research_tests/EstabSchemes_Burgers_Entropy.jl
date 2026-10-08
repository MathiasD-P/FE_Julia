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

function make_parameters(myDG, p, node_name, exact_mass=false)

    param = parameters(
                            pdetype="Burgers",
                            dim=1,
                            bnodes="(" * string(p+1) * ")-"* node_name,
                            qnodes="(" * string(p+1) * ")-"* node_name,
                            enodes="(20)-GLL",
                            fnodes="(1)-GL",
                            refelem="interval",
                            domain="unit_interval_linear",
                            Neldim=20,
                            numfluxtype="LFOLD",
                            ICname="sin_1state",
                            k=2*pi,
                            phi=0.1,
                            av=0.01,
                            BCname="periodic",
                            ODE_solver="SSPRRK33",
                            Nsteps = 20000,
                            calc_entropy=true)
    
    param.dt = 3.0 / param.Nsteps

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

function research_test(p, node_name, DGnames)
    # Construct parameters
    params = Dict()
    for myDG in DGnames
        params[myDG] = make_parameters(myDG, p, node_name)
    end

    paramExact = make_parameters("DGExact", p, node_name)

    # Testing starts here
    outputs = Dict()
    for myDG in DGnames
        println(myDG)
        outputs[myDG] = set_up_and_solve(params[myDG])
        plt = plot(outputs[myDG]["solution"])
        display(plt)
    end

    outputs["DGExact"] = set_up_and_solve(paramExact)

    return outputs
end

function plot_results(outputs, DGnames, ylims)
    Delta = 0.15

    plt_ent = plot(xlabel=L"t",
               ylabel=L"Entropy",
               legend=:topright,
               ylims=ylims,
               dpi = 500)
    
    for myDG in DGnames
        t = collect(range(0, outputs[myDG]["time"], length(outputs[myDG]["entropy"]) + 1))
        t = t[2:end]

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
        end

        ds = sqrt.((t[2:end] .- t[1:end-1]).^2 .+ ((outputs[myDG]["entropy"][2:end] .- outputs[myDG]["entropy"][1:end-1]) ./ 0.25005).^2)
        ds = [0; ds]
        s = [sum(ds[1:i]) for i in 1:length(ds)]
        s[abs.(s) .> 1e10] .= 0
        s[isnan.(s)] .= 0.0
        markerslicer = unique(i -> round.(Int, s ./ Delta)[i], eachindex(s))

        plot!(plt_ent, t, outputs[myDG]["entropy"] ./ 0.25005, label=false, color=mycolor)
        scatter!(plt_ent, t[markerslicer], outputs[myDG]["entropy"][markerslicer] ./ 0.25005, label=mylabel, color=mycolor, markershape=mymarker, markersize=mymarkersize, markerstrokewidth = 0)
    end

    t = collect(range(0, outputs["DGExact"]["time"], length(outputs["DGExact"]["entropy"]) + 1))
    t = t[2:end]

    ds = sqrt.((t[2:end] .- t[1:end-1]).^2 .+ ((outputs["DGExact"]["entropy"][2:end] .- outputs["DGExact"]["entropy"][1:end-1]) ./ 0.25005).^2)
    ds = [0; ds]
    s = [sum(ds[1:i]) for i in 1:length(ds)]
    s[abs.(s) .> 1e10] .= 0
    s[isnan.(s)] .= 0.0
    markerslicer = unique(i -> round.(Int, s ./ Delta)[i], eachindex(s))

    plot!(plt_ent, t, outputs["DGExact"]["entropy"] ./ 0.25005, label=false, color=:orange)
    scatter!(plt_ent, t[markerslicer], outputs["DGExact"]["entropy"][markerslicer] ./ 0.25005, label="Exact Integration", color=:orange, markershape=:rect, markersize=3, markerstrokewidth = 0)

    return plt_ent
end

function plot_solutions(DGnames, t, Nsteps)
    Delta = 0.15

    # Construct parameters
    params = Dict()
    for myDG in DGnames
        params[myDG] = make_parameters(myDG, p, node_name)
        params[myDG].dt = t / Nsteps
        params[myDG].Nsteps = Nsteps
    end

    paramExact = make_parameters("DGExact", p, node_name)
    paramExact.dt = t / Nsteps
    paramExact.Nsteps = Nsteps

    paramRef = make_parameters("DGStd", 1, "GLL")
    paramRef.Neldim=2000
    paramRef.dt = t / 20000
    paramRef.Nsteps = 20000
    paramRef.ODE_solver = "SSPRRK33"

    # Testing starts here
    outputs = Dict()
    for myDG in DGnames
        println(myDG)
        outputs[myDG] = set_up_and_solve(params[myDG])
    end
    outputs["Ref"] = set_up_and_solve(paramRef)
    outputs["DGExact"] = set_up_and_solve(paramExact)

    # Now we plot
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

    chie, _, epts = FE_Julia.extract_volume_quadrature(outputs["DGExact"]["dg"].refelem.bnodestype, make_nodes(paramRef.enodes))
    epts = reduce(vcat, FE_Julia.mapping(outputs["DGExact"]["dg"].mesh, epts, ielem) for ielem in 1:outputs["DGExact"]["dg"].mesh.Nel)
    solution = FE_Julia.block_matmul(chie, outputs["DGExact"]["solution"], outputs["DGExact"]["dg"].mesh.Nel)

    ds = sqrt.((epts[2:end] .- epts[1:end-1]).^2 .+ (solution[2:end] .- solution[1:end-1]).^2)
    ds = [0; ds]
    s = [sum(ds[1:i]) for i in 1:length(ds)]
    s[abs.(s) .> 1e10] .= 0
    markerslicer = unique(i -> round.(Int, s ./ Delta)[i], eachindex(s))

    plot!(plt_sol, epts, solution, label=false, color=:orange)
    scatter!(plt_sol, epts[markerslicer], solution[markerslicer], label="Exact Integration", color=:orange, markershape=:rect, markersize=3, markerstrokewidth = 0)
     
    # We also plot the Ref solution
    chie, _, epts = FE_Julia.extract_volume_quadrature(outputs["Ref"]["dg"].refelem.bnodestype, make_nodes(paramRef.enodes))
    epts = reduce(vcat, FE_Julia.mapping(outputs["Ref"]["dg"].mesh, epts, ielem) for ielem in 1:outputs["Ref"]["dg"].mesh.Nel)

    solution = FE_Julia.block_matmul(chie, outputs["Ref"]["solution"], outputs["Ref"]["dg"].mesh.Nel)
    plot!(plt_sol, epts, solution, label="Reference", color=:black, linestyle=:dot)


    return plt_sol
end

DGnames = ["DGStd", "DGFluxDiff", "DGArtVisc", "DGAddRes"]
p = 5
node_name = "GLL"
outputs = research_test(p, node_name, DGnames)

plt_ent = plot_results(outputs, DGnames, (0,1.1))
display(plt_ent)

plt_sol = plot_solutions(["DGFluxDiff", "DGArtVisc", "DGAddRes"], 3.0, 20000)
display(plt_sol)
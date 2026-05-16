using FE_Julia
using LinearAlgebra
using LaTeXStrings
using Plots

# Default plotting parameters
default(
    fontfamily="Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(12, "Computer Modern"),
    tickfont = font(10, "Computer Modern"),
    legendfont = font(10, "Computer Modern")
)

function make_parameters(bnodes, qnodes, dgtype)
    params = parameters(
                        pdetype="Burgers",
                        dgtype=dgtype,
                        dim=1,
                        bnodes=bnodes,
                        qnodes=qnodes,
                        enodes="(50)-GLL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="LF",
                        twoptfluxtype="EC_split",
                        Rescorr="RescorrEC",
                        AVcoeff="AVdissip",
                        ICname="Burgulence",
                        BCname="periodic",
                        kmax=100)

    return params
end

function norm_Legendre_Vandermonde(x::Vector{Float64}, p::Int64) # evaluate Legendre polynomials on the unit interval [0, 1]
    psi = [x -> 1.0,
           x -> 2.0 * x * sqrt(3),
           x -> 0.5 * (12.0 * x^2 - 1.0) * sqrt(5),
           x -> 0.5 * (40 * x^3 .- 6.0 * x) * sqrt(7),
           x -> (1/8) * (35 * (2.0 * x)^4 - 30 * (2.0 * x)^2 + 3.0) * sqrt(9),
           x -> (1/8) * (63 * (2 * x)^5 - 70 * (2 * x)^3 + 15 * (2 * x)) * sqrt(11)]
    
    V = zeros(length(x), p+1)

    for i in 0:p
        V[:,i+1] .= (psi[i+1]).(x .- 0.5)
    end

    return V
end

function research_test(bnodes, qnodes, qnodesexact, filenames=nothing)
    DGnames = ["DGStd", "DGFluxDiff", "DGArtVisc", "DGAddRes"]

    # We first deal with the exact integration case
    params = make_parameters(bnodes, qnodesexact, "DGStd")

    u0, BChandler, dg = set_up_problem(params)
    residual = similar(u0)
    build_residual!(residual, u0, 0.0, BChandler, dg, params)

    V = norm_Legendre_Vandermonde(dg.refelem.bnodes[:], dg.refelem.Nbnodes - 1)
    exact =  FE_Julia.block_matmul(inv(V), residual, dg.mesh.Nel)

    psi = zeros(dg.refelem.Nbnodes, 4)

    # Now we deal with ES schemes
    for (idg, myDG) in enumerate(DGnames)
        params = make_parameters(bnodes, qnodes, myDG)

        _, BChandler, dg = set_up_problem(params)
        residual = similar(u0)
        build_residual!(residual, u0, 0.0, BChandler, dg, params)

        modalresidual = FE_Julia.block_matmul(inv(V), residual, dg.mesh.Nel)
        psi[:,idg] .= sum(reshape(abs.(modalresidual .- exact), (dg.refelem.Nbnodes, dg.mesh.Nel)), dims=2)[:]
    end

    # we plot state and residual modal decomposition
    enodes = FE_Julia.evaluate(FE_Julia.make_nodes(params.enodes))
    chiinterp = FE_Julia.Lagrange_Vandermonde1D(enodes[:], dg.refelem.bnodes[:])
    epts = reduce(vcat, FE_Julia.mapping(dg.mesh, enodes, ielem) for ielem in 1:dg.mesh.Nel)

    plt = plot(epts, FE_Julia.block_matmul(chiinterp, u0, dg.mesh.Nel), legend=nothing, xlabel=L"x", ylabel=L"u", xlims=(-0.5,0.5), color=:red, dpi=500)

    if isnothing(filenames)
        display(plt)
    else
        savefig(plt, filenames[1] * ".png")
    end

    barnames = [LaTeXString("\\psi_{" * string(i) * "}") for i in 0:(dg.DOF - 1)]
    barnum = collect(1:5)
    w = 0.15
    plt=bar(barnum .- 1.5w, psi[:,1], bar_width=w, label="DG", color="#F0E442", legend=:topleft, xlabel="Normalized Legendre mode", ylabel="Average Residual Error", title="Basis: " * bnodes * "; Quad: " * qnodes, dpi=500)
    bar!(barnum .- 0.5w, psi[:,2], bar_width=w, label="Flux Differencing", color="#E69F00")
    bar!(barnum .+ 0.5w, psi[:,3], bar_width=w, label="Artifical Viscosity", color="#56B4E9")
    bar!(barnum .+ 1.5w, psi[:,4], bar_width=w, label="Residual Correction", color="#009E73")
    xticks!(plt, barnum, barnames)

    if isnothing(filenames)
        display(plt)
    else
        savefig(plt, filenames[2] * ".png")
    end

    return psi
end



bnodes = "(5)-GLL"
qnodes = "(5)-GLL"
qnodesexact = "(8)-GL"

research_test(bnodes, qnodes, qnodesexact, ("outputs/bar_GLLstate", "outputs/bar_GLLmodes"))

bnodes = "(5)-GL"
qnodes = "(5)-GL"
qnodesexact = "(8)-GL"

research_test(bnodes, qnodes, qnodesexact, ("outputs/bar_GLstate", "outputs/bar_GLmodes"))
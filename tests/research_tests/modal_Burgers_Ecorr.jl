using FE_Julia
using LinearAlgebra
using LaTeXStrings
using Plots

# Default plotting parameters
default(
    fontfamily="Computer Modern",
    titlefont = font(10, "Computer Modern"),
    guidefont = font(10, "Computer Modern"),
    tickfont = font(10, "Computer Modern"),
    legendfont = font(8, "Computer Modern")
)

function make_parameters(bnodes, qnodes)
    params = parameters(
                        pdetype="Burgers",
                        dim=1,
                        bnodes=bnodes,
                        qnodes=qnodes,
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=1,
                        numfluxtype="EC_split",
                        twoptfluxtype="EC_split",
                        Rescorr="RescorrEC",
                        AVcoeff="AVEC",
                        ICname="GassnerBurgers",
                        sourcename="GassnerBurgers",
                        BCname="periodic")

    return params
end

function norm_Legendre_Vandermonde(x::Vector{Float64}, p::Int64) # evaluate Legendre polynomials on the unit interval
    psi = [x -> 1.0,
           x -> 2.0 * x * sqrt(3),
           x -> 0.5 * (12.0 * x^2 - 1.0) * sqrt(5),
           x -> 0.5 * (40 * x^3 .- 6.0 * x) * sqrt(7),
           x -> (1/8) * (35 * (2.0 * x)^4 - 30 * (2.0 * x)^2 + 3.0) * sqrt(9),
           x -> (1/8) * (63 * (2 * x)^5 - 70 * (2 * x)^3 + 15 * (2 * x)) * sqrt(11)]
    
    V = zeros(length(x), p+1)

    for i in 0:p
        V[:,i+1] .= (psi[i+1]).(x)
    end

    return V
end

function research_test(bnodes, qnodes, mydg)
    params = make_parameters(bnodes, qnodes)

    # Standard DG residual
    params.dgtype = "DGArtVisc"
    params.AVcoeff = "NoAV"
    u0, BChandler, dg = set_up_problem(params)
    stdresidual = similar(u0)
    build_residual!(stdresidual, u0, 0.0, BChandler, dg, params)
    params.AVcoeff = "AVEC"

    # ES DG residual
    params.dgtype = mydg
    u0, BChandler, dg = set_up_problem(params)
    residual = similar(u0)
    build_residual!(residual, u0, 0.0, BChandler, dg, params)

    dres = norm_Legendre_Vandermonde(dg.bpts[:], dg.DOF - 1) \ (residual .- stdresidual)

    # Now, we plot
    myplots = []
    barnames = [LaTeXString("\\phi_{" * string(i) * "}") for i in 0:(dg.DOF - 1)]
    for s in 1:dg.Nstates
        plt=plot(bar(barnames, dres[:,s] ./ norm(dres[:,s])), ylabel="Normalized Magnitude", ylims=(-1,1), label=nothing, legend=false)
        push!(myplots, plt)
    end
    plt = plot(myplots..., layout=(1,dg.Nstates), plot_title= params.pdetype * "; Basis: " * bnodes * "; Quad: " * qnodes)
    display(plt)

    return dres
end

bnodes = "(6)-GLL"
qnodes = "(6)-GLL"
mydg = "DGAddRes"

dres = research_test(bnodes, qnodes, mydg)

bnodes = "(6)-GLL"
qnodes = "(6)-GL"
mydg = "DGAddRes"

dres = research_test(bnodes, qnodes, mydg)

bnodes = "(6)-GLL"
qnodes = "(7)-GL"
mydg = "DGAddRes"

dres = research_test(bnodes, qnodes, mydg)
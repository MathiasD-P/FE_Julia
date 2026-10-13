include("linearized_VonNeumann.jl")

using Plots
using LaTeXStrings

default(
    fontfamily="Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(12, "Computer Modern"),
    tickfont = font(10, "Computer Modern"),
    legendfont = font(10, "Computer Modern")
)

function research_test(prange, ref, physics, ic, nodetype, vglims)
    plt = plot(xlabel = L"\theta / (p+1)", ylabel=L"Re(\omega) / (p+1)", title = nodetype * " nodes", dpi=700)
    plt1 = plot(xlabel = L"\theta / (p+1)", ylabel=L"Im(\omega) / (p+1)", title = nodetype * " nodes", dpi=700)
    plt2 = plot(xlabel = L"\theta / (p+1)", ylabel=L"v_g", title = nodetype * " nodes", xlims=(0,1), ylims=vglims, dpi=700)

    for p in prange
        qnodes = FE_Julia.Tensorprod_nodes(nodetype, p+1, 1)
        bnodes = FE_Julia.Tensorprod_nodes("GLL", p+1, 1)

        dgstd = DGStd(1, RefElemStd(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
        dgfdf = DGFluxDiff(1, RefElemSBP(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
        dgavi = DGArtVisc(1, RefElemStd(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
        dgadr = DGAddRes(1, RefElemStd(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))

        outstd = linearized_VN(dgstd, physics, ic.(dgstd.bpts))
        outfdf = linearized_VN(dgfdf, physics, ic.(dgfdf.bpts))
        outavi = linearized_VN(dgavi, physics, ic.(dgavi.bpts))
        outadr = linearized_VN(dgadr, physics, ic.(dgadr.bpts))

        plot!(plt, outstd[1] ./ (p+1), real(outstd[2])./ (p+1), color=:black, label=false)
        plot!(plt, outfdf[1] ./ (p+1), real(outfdf[2])./ (p+1), color=:red, label=false)
        plot!(plt, outavi[1] ./ (p+1), real(outavi[2])./ (p+1), color=:blue, label=false)
        plot!(plt, outadr[1] ./ (p+1), real(outadr[2])./ (p+1), color=:green, label=false)

        plot!(plt1, outstd[1] ./ (p+1), imag(outstd[2])./ (p+1), color=:blue)
        plot!(plt1, outfdf[1] ./ (p+1), imag(outfdf[2])./ (p+1), color=:red)

        wavenum, groupvel = compute_groupvel(outstd[1], real(outstd[2]))
        plot!(plt2, wavenum./ (p+1), groupvel, color=:black, label=false)
        wavenum, groupvel = compute_groupvel(outfdf[1], real(outfdf[2]))
        plot!(plt2, wavenum./ (p+1), groupvel, color=:red, label=false)
        wavenum, groupvel = compute_groupvel(outavi[1], real(outavi[2]))
        plot!(plt2, wavenum./ (p+1), groupvel, color=:blue, label=false)
        wavenum, groupvel = compute_groupvel(outadr[1], real(outadr[2]))
        plot!(plt2, wavenum./ (p+1), groupvel, color=:green, label=false)
    end

    plot!(plt, [0,2], [0, 2*ref], color=:black, linestyle=:dot, label=false)
    plot!(plt2, [0,1], [ref, ref], color=:black, linestyle=:dot, label=false)

    return plt, plt2
end

# constant state

# physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.0], -2=>[1.0]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
# ic = x -> 1.0
# ref = 1.0
# nodetype = "GLL"
# prange = 4:4
# vglims = (0.99, 1)
# pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

# display(pltdisp)
# display(pltvg)

# physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.0], -2=>[1.0]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
# ic = x -> 1.0
# ref = 1.0
# nodetype = "GL"
# prange = 4:4
# vglims = (0.99, 1)
# pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

# display(pltdisp)
# display(pltvg)


# linear state

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 0.910239226627
nodetype = "GLL"
prange = 3:5
vglims = (0.88,0.95)
pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

display(pltdisp)
display(pltvg)

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 0.910239226627
nodetype = "GL"
prange = 3:5
vglims = (0.88,0.95)
pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

display(pltdisp)
display(pltvg)


# # Cubic state

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = 0.960116273343
nodetype = "GLL"
prange = 4:4
vglims = (0.92, 1.01)
pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

display(pltdisp)
display(pltvg)

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = 0.960116273343
nodetype = "GL"
prange = 4:4
vglims = (0.92, 1.01)
pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

display(pltdisp)
display(pltvg)


# # sin state

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = 0.866025403784
nodetype = "GLL"
prange = 4:4
vglims = (0.84, 0.88)
pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

display(pltdisp)
display(pltvg)


physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = 0.866025403784
nodetype = "GL"
prange = 4:4
vglims = (0.84, 0.88)
pltdisp, pltvg = research_test(prange, ref, physics, ic, nodetype, vglims)

display(pltdisp)
display(pltvg)

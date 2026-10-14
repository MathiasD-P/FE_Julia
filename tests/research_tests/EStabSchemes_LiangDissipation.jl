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

function research_test(prange, ref, physics, ic, nodetype)
    plt = plot(xlabel = L"\theta / (p+1)", ylabel=L"\omega^* - \omega_c", title = nodetype * " nodes", dpi=700)

    physicsref = deepcopy(physics)
    physicsref.BChandler = Dict(-1=> [1.0], -2=>[1.0])
    icref = x -> 1.0

    for p in prange
        qnodes = FE_Julia.Tensorprod_nodes(nodetype, p+1, 1)
        bnodes = FE_Julia.Tensorprod_nodes("GLL", p+1, 1)

        dgstd = DGStd(1, RefElemStd(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
        dgfdf = DGFluxDiff(1, RefElemSBP(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
        dgavi = DGArtVisc(1, RefElemStd(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
        dgadr = DGAddRes(1, RefElemStd(bnodes, qnodes, make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))

        outstd = linearized_nonmodal_LIANG(dgstd, physics, ic.(dgstd.bpts))
        outstdref = linearized_nonmodal(dgstd, physicsref, icref.(dgstd.bpts))
        outfdf = linearized_nonmodal(dgfdf, physics, ic.(dgfdf.bpts))
        outfdfref = linearized_nonmodal(dgfdf, physicsref, icref.(dgfdf.bpts))
        outavi = linearized_nonmodal(dgavi, physics, ic.(dgavi.bpts))
        outaviref = linearized_nonmodal(dgavi, physicsref, icref.(dgavi.bpts))
        outadr = linearized_nonmodal(dgadr, physics, ic.(dgadr.bpts))
        outadrref = linearized_nonmodal(dgadr, physicsref, icref.(dgadr.bpts))

        plot!(plt, outstd[1][2:end] ./ (p+1), real(outstd[2] - outstdref[2])[2:end], color=:black, label="DG " * L"(p="*string(p) *")")
        plot!(plt, outfdf[1][2:end] ./ (p+1), real(outfdf[2] - outfdfref[2])[2:end], color=:red, label="FD " * L"(p="*string(p) *")", linestyle=:dashdot)
        plot!(plt, outavi[1][2:end] ./ (p+1), real(outavi[2] - outaviref[2])[2:end], color=:blue, label=label="AV " * L"(p="*string(p) *")", linestyle=:dashdotdot)
        plot!(plt, outadr[1][2:end] ./ (p+1), real(outadr[2] - outadrref[2])[2:end], color=:green, label=label="RC " * L"(p="*string(p) *")", linestyle=:dash)
    end

    plot!(plt, [0,pi], [ref, ref], color=:black, linestyle=:dot, label=false)

    return plt
end

# Linear state

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 1.0
nodetype = "GLL"
prange = 3:3
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/linear_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 1.0
nodetype = "GL"
prange = 3:3
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/linear_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 1.0
nodetype = "GLL"
prange = 4:4
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/linear_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 1.0
nodetype = "GL"
prange = 4:4
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/linear_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 1.0
nodetype = "GLL"
prange = 5:5
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/linear_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [1.5], -2=>[0.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 - 0.5 * (2*x-1)
ref = 1.0
nodetype = "GL"
prange = 5:5
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/linear_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")


# Cubic state

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = -1.0
nodetype = "GLL"
prange = 3:3
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/cubic_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = -1.0
nodetype = "GL"
prange = 3:3
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/cubic_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = -1.0
nodetype = "GLL"
prange = 4:4
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/cubic_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = -1.0
nodetype = "GL"
prange = 4:4
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/cubic_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = -1.0
nodetype = "GLL"
prange = 5:5
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/cubic_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5 * (2*x-1)^3
ref = -1.0
nodetype = "GL"
prange = 5:5
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/cubic_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

# sin state

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = -1.0
nodetype = "GLL"
prange = 3:3
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/sin_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = -1.0
nodetype = "GL"
prange = 3:3
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/sin_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = -1.0
nodetype = "GLL"
prange = 4:4
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/sin_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = -1.0
nodetype = "GL"
prange = 4:4
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/sin_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = -1.0
nodetype = "GLL"
prange = 5:5
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/sin_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")

physics = FE_Julia.PhysProp(FE_Julia.Burgers(1), nothing, Dict(-1=> [0.5], -2=>[1.5]), FE_Julia.LFNumFlux(), FE_Julia.ECSplitTPFlux(), FE_Julia.AVdissip(nothing), FE_Julia.ResCorrEC())
ic = x -> 1.0 + 0.5*sin(0.5*pi*(2*x-1))
ref = -1.0
nodetype = "GL"
prange = 5:5
pltdissip = research_test(prange, ref, physics, ic, nodetype)

display(pltdissip)

savefig(pltdissip, "outputs/sin_p" * string(prange[1]) * "_dissip_" *nodetype * ".pdf")
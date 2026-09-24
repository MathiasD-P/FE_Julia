using FE_Julia
using LinearAlgebra
using Plots
using LaTeXStrings

default(
    fontfamily="Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(12, "Computer Modern"),
    tickfont = font(10, "Computer Modern"),
    legendfont = font(10, "Computer Modern")
)

"""
    classical_VN(p::Int64, numflux::NumFlux, theta::Vector{Float64})

Performs classical Von Neumann analysis for standard DG.
"""
function classical_VN(refelem::RefElemStd, numflux::NumFlux, delta=0.001, all::Bool=false)
    PDE = FE_Julia.LinAdv(1, 1.0)
    dg = DGStd(1, refelem, FE_Julia.make_interval([0, 1.0], [-1, -2]))

    if all
        wavenum = collect(0.0:delta:(refelem.Nbnodes*pi))
        omega = zeros(ComplexF64, length(wavenum), refelem.Nbnodes)
        modes = zeros(ComplexF64, length(omega), refelem.Nbnodes, refelem.Nbnodes)
    else
        wavenum = collect(0.0:delta:(refelem.Nbnodes*pi))
        omega = zeros(ComplexF64, size(wavenum))
        modes = zeros(ComplexF64, length(omega), refelem.Nbnodes)
    end

    domega = 1.0

    for thetai in 2:length(wavenum)
        theta = wavenum[thetai]
        ext_mat = [refelem.chif[2,:]' .* exp(-im * theta); refelem.chif[1,:]' .* exp(-im * theta)]
        VN_mat = -refelem.Dh[1] .+ refelem.LIFT[1] * (refelem.chif .- FE_Julia.compute_numflux(refelem.chif, ext_mat, repeat(dg.nphys, 1, refelem.Nbnodes), numflux, PDE)[1]) # haha, weird overloading of numflux

        # Solve eigenvalue problem
        F = eigen(VN_mat)

        if all
            omega[thetai, :] .= F.values
            modes[thetai,:,:] .= F.vectors
        else
            # Find the physical mode
            targeti = argmin(abs.(F.values .- (omega[thetai-1] + domega * delta)))
            omega[thetai] = F.values[targeti]
            modes[thetai,:] = F.vectors[:,targeti] .* (conj(F.vectors[1,targeti]) / abs(F.vectors[1,targeti])) # we make sure to align phase as well

            # compute linear approximation of omega
            domega = (omega[thetai]-omega[thetai-1]) / delta
        end
    end

    omega .= omega ./ -im

    return wavenum, omega, modes
end


function L2_wave(dg::FE_Julia.DG, physics::FE_Julia.PhysProp, wavenum::Vector{Float64}, mean::Float64, NQ::Int64=1000)
    # Overwrite mesh in DG object
    dg = typeof(dg)(1, dg.refelem, FE_Julia.make_interval([0, 1.0], [-1, -2]))

    # Construct refelem for L2 projection
    refelem_proj = RefElemStd(dg.refelem.bnodestype, FE_Julia.Tensorprod_nodes("GL", NQ, 1), dg.refelem.fnodestype)

    mu = similar(wavenum) # dissipation coefficient
    nu = similar(wavenum) # dispersion coefficient

    for itheta in eachindex(wavenum)
        theta = wavenum[itheta]

        # Start by constructing the relevant L2 projections
        sinL2 = refelem_proj.Ph * sin.(theta .* refelem_proj.qnodes) .+ mean
        cosL2 = refelem_proj.Ph * cos.(theta .* refelem_proj.qnodes) .+ mean

        # Construct the associated BCHandler
        sinBC = Dict(-1 => [sin(0.0) + mean], -2 => [sin(theta) + mean])
        cosBC = Dict(-1 => [cos(0.0) + mean], -2 => [cos(theta) + mean])

        # construct specific physics
        physics_sin = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, sinBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)
        physics_cos = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, cosBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)

        # now we construct residuals
        res_sin = similar(sinL2)
        build_residual!(res_sin, sinL2, 0.0, dg, physics_sin)
        res_cos = similar(cosL2)
        build_residual!(res_cos, cosL2, 0.0, dg, physics_cos)

        # we attempt to extract dispersion and dissipation
        sinL2 .= sinL2 .- mean
        cosL2 .= cosL2 .- mean
        mu[itheta] = sum((res_sin .* sinL2 .+ res_cos .* cosL2) ./ (sinL2.^2 .+ cosL2.^2)) / dg.refelem.Nbnodes
        nu[itheta] = -sum((res_sin .* cosL2 .- res_cos .* sinL2) ./ (sinL2.^2 .+ cosL2.^2)) / dg.refelem.Nbnodes
    end

    return mu, nu
end


function modified_VN(dg::FE_Julia.DG, physics::FE_Julia.PhysProp, mean::Float64, delta::Float64=0.001)
    # Overwrite mesh in DG object
    dg = typeof(dg)(1, dg.refelem, FE_Julia.make_interval([0, 1.0], [-1, -2]))

    # Compute VN modes for standard (linear) DG
    wavenum, omega, modes = classical_VN(RefElemStd(dg.refelem.bnodestype, dg.refelem.qnodestype, dg.refelem.fnodestype), physics.numflux, delta)

    # Now we substitute in the residuals
    gamma = similar(omega) # dispersion-dissipation relation

    for thetai in eachindex(wavenum)
        theta = wavenum[thetai]

        a = reshape(real(modes[thetai,:]), dg.refelem.Nbnodes, 1)
        b = reshape(imag(modes[thetai,:]), dg.refelem.Nbnodes, 1)

        a0 = (dg.refelem.chif * a)[1]
        a1 = (dg.refelem.chif * a)[2]
        b0 = (dg.refelem.chif * b)[1]
        b1 = (dg.refelem.chif * b)[2]

        u = a .+ mean
        w = b .+ mean

        uBC = Dict(-1 => [a1 * cos(theta) + b1 * sin(theta) + mean], -2 => [a0 * cos(theta) - b0 * sin(theta) + mean])
        wBC = Dict(-1 => [-a1 * sin(theta) + b1 * cos(theta) + mean], -2 => [a0 * sin(theta) + b0 * cos(theta) + mean])

        uphysics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, uBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)
        wphysics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, wBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)

        ures = similar(u)
        build_residual!(ures, u, 0.0, dg, uphysics)
        wres = similar(w)
        build_residual!(wres, w, 0.0, dg, wphysics)

        # Now we assemble and compute gamma
        v = ures .+ im .* wres
        gamma[thetai] = sum(v ./ modes[thetai,:]) / length(v) / -im # here we compute gamma using the average
    end

    return wavenum, gamma
end

function modified_VN_spurious(dg::FE_Julia.DG, physics::FE_Julia.PhysProp, mean::Float64, delta::Float64=0.001)
    # Overwrite mesh in DG object
    dg = typeof(dg)(1, dg.refelem, FE_Julia.make_interval([0, 1.0], [-1, -2]))

    # Compute VN modes for standard (linear) DG
    wavenum, omega, modes = classical_VN(RefElemStd(dg.refelem.bnodestype, dg.refelem.qnodestype, dg.refelem.fnodestype), physics.numflux, delta)

    # Now we substitute in the residuals
    gamma = similar(omega) # dispersion-dissipation relation
    kappa = similar(omega) # spurious mode activation

    for thetai in eachindex(wavenum)
        theta = wavenum[thetai]

        a = reshape(real(modes[thetai,:]), dg.refelem.Nbnodes, 1)
        b = reshape(imag(modes[thetai,:]), dg.refelem.Nbnodes, 1)

        a0 = (dg.refelem.chif * a)[1]
        a1 = (dg.refelem.chif * a)[2]
        b0 = (dg.refelem.chif * b)[1]
        b1 = (dg.refelem.chif * b)[2]

        u = a .+ mean
        w = b .+ mean

        uBC = Dict(-1 => [a1 * cos(theta) + b1 * sin(theta) + mean], -2 => [a0 * cos(theta) - b0 * sin(theta) + mean])
        wBC = Dict(-1 => [-a1 * sin(theta) + b1 * cos(theta) + mean], -2 => [a0 * sin(theta) + b0 * cos(theta) + mean])

        uphysics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, uBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)
        wphysics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, wBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)

        if dg isa DGArtVisc # Reconstruct physics object with our gradient variables.
            ux = FE_Julia.build_sigma(u, 0.0, dg, uphysics)[1]
            wx = FE_Julia.build_sigma(w, 0.0, dg, uphysics)[1]

            ux0 = (dg.refelem.chif * ux)[1]
            ux1 = (dg.refelem.chif * ux)[2]
            wx0 = (dg.refelem.chif * wx)[1]
            wx1 = (dg.refelem.chif * wx)[2]

            uxBC = (Dict(-1 => [ux1 * cos(theta) + wx1 * sin(theta)], -2 => [ux0 * cos(theta) - wx0 * sin(theta)]),)
            wxBC = (Dict(-1 => [-ux1 * sin(theta) + wx1 * cos(theta)], -2 => [ux0 * sin(theta) + wx0 * cos(theta)]),)

            uphysics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, uBC, uxBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)
            wphysics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, wBC, wxBC, physics.numflux, physics.tpflux, physics.artvisc, physics.rescorr)
        end

        ures = similar(u)
        build_residual!(ures, u, 0.0, dg, uphysics)
        wres = similar(w)
        build_residual!(wres, w, 0.0, dg, wphysics)

        # Now we assemble and compute gamma
        v = vec(ures) .+ im .* vec(wres)
        gamma[thetai] = (dot(modes[thetai,:], v) / norm(modes[thetai,:])^2) / -im

        kappa[thetai] = norm(v .- (dot(modes[thetai,:], v) / norm(modes[thetai,:])^2) .* modes[thetai,:]) / norm(v)
    end

    return wavenum, gamma, kappa
end

### TESTS FOR CLASSICAL VON NEUMANN

# refe = RefElemStd(make_nodes("(4)-GL"), make_nodes("(4)-GL"), make_nodes("interval", "(1)-GLL"))
# numflux = FE_Julia.UpwindNumFlux()

# out = classical_VN(refe, numflux, 0.01)
# plt1=plot(out[1], imag(out[2]))
# display(plt1)
# plt2=plot(out[1], real(out[2]))
# display(plt2)


# ### TESTS FOR L2_WAVE

# dg = DGStd(1, RefElemStd(make_nodes("(4)-GLL"), make_nodes("(4)-GL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
# physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), nothing, nothing, nothing)
# wavenum = collect(0.0:0.01:(4*pi))
# mean = 20.0

# mu, nu = L2_wave(dg, physics, wavenum, mean)

# plt1 = plot(wavenum,mu)
# plt2 = plot(wavenum,nu)


# dg = DGFluxDiff(1, RefElemSBP(make_nodes("(4)-GLL"), make_nodes("(4)-GL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
# physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), FE_Julia.LogMeanTPFlux(), nothing, nothing)
# wavenum = collect(0.0:0.01:(4*pi))
# mean = 1.5

# mu, nu = L2_wave(dg, physics, wavenum, mean)

# plot!(plt1,wavenum,mu)
# plot!(plt2,wavenum,nu)


# dg = DGAddRes(1, RefElemStd(make_nodes("(4)-GLL"), make_nodes("(4)-GL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
# physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), nothing, nothing, FE_Julia.ResCorrEC())
# wavenum = collect(0.0:0.01:(4*pi))
# mean = 1.5

# mu, nu = L2_wave(dg, physics, wavenum, mean)

# println(mu)

# plot!(plt1,wavenum,mu)
# plot!(plt2,wavenum,nu)
# plot!(plt1, wavenum, zeros(size(wavenum)), color=:black)
# plot!(plt2, wavenum, wavenum, color=:black)
# display(plt1)
# display(plt2)


### TESTS FOR MODIFED VON NEUMANN

# mean = 0.5

# dg = DGStd(1, RefElemStd(make_nodes("(4)-GLL"), make_nodes("(4)-GL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
# physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), nothing, nothing, nothing)

# wavenum, gamma = modified_VN(dg, physics, mean)

# plt1 = plot(wavenum, real(gamma)) # dispersion
# plt2 = plot(wavenum, imag(gamma)) # dissipation


# dg = DGFluxDiff(1, RefElemSBP(make_nodes("(4)-GLL"), make_nodes("(4)-GL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
# physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), FE_Julia.LogMeanTPFlux(), nothing, nothing)

# wavenum, gamma = modified_VN(dg, physics, mean)

# plot!(plt1, wavenum, real(gamma)) # dispersion
# plot!(plt2, wavenum, imag(gamma)) # dissipation

# dg = DGAddRes(1, RefElemStd(make_nodes("(4)-GLL"), make_nodes("(5)-GL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
# physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), nothing, nothing, FE_Julia.ResCorrDissip())

# wavenum, gamma = modified_VN(dg, physics, mean)

# plot!(plt1, wavenum, real(gamma)) # dispersion
# plot!(plt2, wavenum, imag(gamma)) # dissipation

# display(plt1)
# display(plt2)

# ### TESTS FOR MODIFED VON NEUMANN with spurious

mean = 1.0

dg = DGStd(1, RefElemStd(make_nodes("(5)-GLL"), make_nodes("(5)-GLL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), nothing, nothing, nothing)

wavenum, gamma, kappa = modified_VN_spurious(dg, physics, mean)

plt1 = plot(wavenum, real(gamma), xlabel=L"\theta", ylabel=L"Re(\omega)", label="Standard", color=:black) # dispersion
plt2 = plot(wavenum, imag(gamma), xlabel=L"\theta", ylabel=L"Im(\omega)", label="Standard", color=:black) # dissipation
plt3 = plot(wavenum, real(kappa), xlabel=L"\theta", ylabel=L"\kappa", label="Standard", color=:black)


dg = DGFluxDiff(1, RefElemSBP(make_nodes("(5)-GLL"), make_nodes("(5)-GLL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), FE_Julia.LogMeanTPFlux(), nothing, nothing)

wavenum, gamma, kappa = modified_VN_spurious(dg, physics, mean)

plot!(plt1, wavenum, real(gamma), label="Flux Diff.", color=:red) # dispersion
plot!(plt2, wavenum, imag(gamma), label="Flux Diff.", color=:red) # dissipation
plot!(plt3, wavenum, real(kappa), label="Flux Diff.", color=:red)

dg = DGAddRes(1, RefElemStd(make_nodes("(5)-GLL"), make_nodes("(5)-GLL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), nothing, nothing, FE_Julia.ResCorrEC())

wavenum, gamma, kappa = modified_VN_spurious(dg, physics, mean)

plot!(plt1, wavenum, real(gamma), label="Res. Corr.", color=:green) # dispersion
plot!(plt2, wavenum, imag(gamma), label="Res. Corr.", color=:green) # dissipation
plot!(plt3, wavenum, real(kappa), label="Res. Corr.", color=:green)

dg = DGArtVisc(1, RefElemStd(make_nodes("(5)-GLL"), make_nodes("(5)-GLL"), make_nodes("interval", "(1)-GLL")), FE_Julia.make_interval([0, 1.0], [-1, -2]))
physics = FE_Julia.PhysProp(FE_Julia.LinAdvLogE(1, 1.0), nothing, Dict(), FE_Julia.UpwindNumFlux(), nothing, FE_Julia.AVdissip(nothing), nothing)

wavenum, gamma, kappa = modified_VN_spurious(dg, physics, mean)

plot!(plt1, wavenum, real(gamma), label="Art. Visc.", color=:blue) # dispersion
plot!(plt2, wavenum, imag(gamma), label="Art. Visc.", color=:blue) # dissipation
plot!(plt3, wavenum, real(kappa), label="Art. Visc.", color=:blue)

display(plt1)
display(plt2)
display(plt3)
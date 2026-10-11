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

function linearized_VN(dg::Union{DGStd, DGFluxDiff, DGAddRes, DGArtVisc}, physics::FE_Julia.PhysProp, u0::Matrix{Float64}, delta=0.001)
    # First, we rebuild a DG object on the reference element
    dg = typeof(dg)(1, dg.refelem, FE_Julia.make_interval([0, 1.0], [-1, -2]))

    # preallocate measured quantities
    wavenum = collect(0.0:delta:(dg.refelem.Nbnodes*pi))
    omega = zeros(ComplexF64, size(wavenum))
    modes = zeros(ComplexF64, length(omega), dg.refelem.Nbnodes)

    domega = -1.0 * im

    for thetai in 2:length(wavenum)
        theta = wavenum[thetai]

        # Assemble VN matrix
        jac = FE_Julia.residual_jacobian_Dirichlet(u0, 0.0, dg, physics)
        VN_mat = reduce(vcat, transpose((jacrow[1] * exp(-im * theta)) .* dg.refelem.chif[2,:] .+ (jacrow[2] * exp(im * theta)) .* dg.refelem.chif[1,:] .+ jacrow[3:end]) for jacrow in eachrow(jac))

        # Solve eigenvalue problem
        F = eigen(VN_mat)

        # Find the physical mode
        targeti = argmin(abs.(F.values .- (omega[thetai-1] + domega * delta)))
        omega[thetai] = F.values[targeti]
        modes[thetai,:] = F.vectors[:,targeti] .* (conj(F.vectors[1,targeti]) / abs(F.vectors[1,targeti])) # we make sure to align phase as well

        # compute linear approximation of omega
        domega = (omega[thetai]-omega[thetai-1]) / delta
    end

    omega .= omega ./ -im

    return wavenum, omega, modes
end

function compute_groupvel(wavenum, omega) # just central differences on dispersion curve
    groupvel = (omega[3:end] .- omega[1:end-2]) ./ (wavenum[3:end] .- wavenum[1:end-2])

    return wavenum[2:end-1], groupvel
end

function linearized_nonmodal(dg::Union{DGStd, DGFluxDiff, DGAddRes}, physics::FE_Julia.PhysProp, u0::Matrix{Float64}, delta::Float64=0.001, NQ::Int64=1000)
    # First, we rebuild a DG object on the reference element
    dg = typeof(dg)(1, dg.refelem, FE_Julia.make_interval([0, 1.0], [-1, -2]))
    refelem_proj = RefElemStd(dg.refelem.bnodestype, FE_Julia.Tensorprod_nodes("GL", NQ, 1), dg.refelem.fnodestype)

    # preallocate measured quantities
    wavenum = collect(0.0:delta:(dg.refelem.Nbnodes*pi))
    omega = zeros(ComplexF64, size(wavenum))

    for thetai in 2:length(wavenum)
        theta = wavenum[thetai]

        # Assemble VN matrix
        jac = FE_Julia.residual_jacobian_Dirichlet(u0, 0.0, dg, physics)
        VN_mat = reduce(vcat, transpose((jacrow[1] * exp(-im * theta)) .* dg.refelem.chif[2,:] .+ (jacrow[2] * exp(im * theta)) .* dg.refelem.chif[1,:] .+ jacrow[3:end]) for jacrow in eachrow(jac))

       evec = refelem_proj.Ph * exp.(im .* theta .* refelem_proj.qnodes)
       omega[thetai] = real(evev' * VN_mat * evec) / (evec' * dg.relem.M * evec)
    end

    return omega
end
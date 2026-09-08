#####################################################################
# Everything concerning the problem's physics
#####################################################################

# THE SOLUTION IS ALWAYS A MATRIX OF SIZE (Nstate, NDOF)

#####################################################################
# Governing Equations
#####################################################################

# Equations supported by our solvers

abstract type GoverningPDE end

struct LinAdv <: GoverningPDE
    a::Float64
    dim::Int64
    Nstates::Int64

    function LinAdv(dim::Int, a::Real)
        new(
            a,
            dim,
            1
        )
    end
end

struct LinAdvLogE <: GoverningPDE
    a::Float64
    dim::Int64
    Nstates::Int64

    function LinAdvLogE(dim::Int, a::Real)
        new(
            a,
            dim,
            1
        )
    end
end

struct Burgers <: GoverningPDE
    dim::Int64
    Nstates::Int64

    function Burgers(dim::Int)
        new(
            dim,
            1
        )
    end
end

struct EulerPerfGas <: GoverningPDE
    gamma::Float64
    dim::Int64
    Nstates::Int64

    function EulerPerfGas(dim::Int, gamma::Real)
        new(
            gamma,
            dim,
            2 + dim
        )
    end
end

# Physical fluxes for our supported PDEs

function compute_physflux(u::AbstractMatrix, PDE::LinAdv)
    return (PDE.a .* u,)
end

function compute_physflux(u::AbstractMatrix, PDE::LinAdvLogE)
    return (PDE.a .* u,)
end

function compute_physflux(u::AbstractMatrix, PDE::Burgers)
    return (0.5 .* u.^2,)
end

function compute_physflux(u::AbstractMatrix, PDE::EulerPerfGas)
    return Euler_physflux(u, PDE)
end


# Conservative variables to entropy variables for our supported PDEs

function compute_evar(u::AbstractMatrix, PDE::LinAdvLogE)
    return -u.^(-1)
end

function compute_evar(u::AbstractMatrix, PDE::Burgers)
    return u
end

function compute_evar(u::AbstractMatrix, PDE::EulerPerfGas)
    return Euler_evar(u, PDE)
end


# Entropy variables to conservative variables for our supported PDEs

function compute_cvar(v::AbstractMatrix, PDE::LinAdvLogE)
    return -v.^(-1)
end

function compute_cvar(v::AbstractMatrix, PDE::Burgers)
    return v
end

function compute_cvar(v::AbstractMatrix, PDE::EulerPerfGas)
    return Euler_cvar(v, PDE)
end

# Compute the entropy for our supported PDEs

function compute_local_entropy(u::AbstractMatrix, PDE::LinAdv)
    return 0.5 .* u.^2
end

function compute_local_entropy(u::AbstractMatrix, PDE::LinAdvLogE)
    return -log.(u)
end

function compute_local_entropy(u::AbstractMatrix, PDE::Burgers)
    return 0.5 .* u.^2
end

function compute_local_entropy(u::AbstractMatrix, PDE::EulerPerfGas)
    return -u[:,1] .* Euler_cvar_entropy(u, PDE)
end

# Compute the entropy potential for our supported PDEs

function compute_cvar_potential(u::AbstractArray, PDE::LinAdvLogE)
    return (PDE.a .* log.(u),)
end

function compute_cvar_potential(u::AbstractArray, PDE::Burgers)
    return ((1/6) .* u.^3,)
end

function compute_cvar_potential(u::AbstractArray, PDE::EulerPerfGas)
    return Tuple((PDE.gamma-1) .* u[:,istate] for istate in 2:1+PDE.dim) # CAREFUL, MISTAKE IN (Chan, 2025)
end


# Compute Hessian from conservative variables
# Operates on one node at the time.

function compute_cvar_Hessian(u::AbstractVector, PDE::Burgers)
    K = Matrix{Float64}(undef, 1, 1)
    K[1,1] = 1.0
    return K
end

function compute_cvar_Hessian(u::AbstractVector, PDE::EulerPerfGas)
    return Euler_cvar_Hessian(u, PDE)
end


#####################################################################
# Source Terms
#####################################################################

# The specific source expressions are defined in BCsICsSources.jl

abstract type CLawSource end


#####################################################################
# Numerical Fluxes
#####################################################################

abstract type NumFlux end

struct CentralNumFlux <: NumFlux end

struct UpwindNumFlux <: NumFlux end

struct LFNumFlux <: NumFlux end

struct ECSplitNumFlux <: NumFlux end

struct ECChandrashekarNumFlux <: NumFlux end

struct ESChandrashekarDissipNumFlux <: NumFlux end


# Linear advection numerical fluxes

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys, numflux::CentralNumFlux, PDE::Union{LinAdv, LinAdvLogE})
    return (0.5 .* PDE.a .* (up .+ un),)
end

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::AbstractMatrix, numflux::UpwindNumFlux, PDE::Union{LinAdv, LinAdvLogE})
    f = copy(up)

    if PDE.a * nphys[1] >= 0
        index = nphys .== nphys[1]
    else
        index = nphys .!= nphys[1]
    end

    f[index] = un[index]
    return (PDE.a .* f,)
end


# Burgers numerical fluxes

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::Union{AbstractMatrix,Nothing}, numflux::CentralNumFlux, PDE::Burgers)
    return (0.25 .* (up.^2 .+ un.^2),)
end

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::AbstractMatrix, numflux::UpwindNumFlux, PDE::Burgers)
    f = similar(up)

    vel = 0.5 .* (un.^2 .- up.^2) ./ (un .- up) # local velocity
    slicing = vel .* nphys .>= 0
    notslicing = .!slicing
    f[slicing] .= 0.5 .* un[slicing].^2
    f[notslicing] .= 0.5 .* up[notslicing].^2

    return (f,)
end

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::Union{AbstractMatrix,Nothing}, numflux::LFNumFlux, PDE::Burgers)
    return (0.25 .* ((un.^2 .+ up.^2) .- max.(abs.(up), abs.(un)) .* (up .- un) .* nphys),)
end

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::Union{AbstractMatrix,Nothing}, numflux::ECSplitNumFlux, PDE::Burgers)
    return ((1/6) .* (un.^2 .+ up .* un .+ up.^2),)
end


# Euler numerical fluxes

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::Union{AbstractMatrix,Nothing}, numflux::CentralNumFlux, PDE::EulerPerfGas)
    fp = Euler_physflux(up, PDE)
    fn = Euler_physflux(up, PDE)

    for dir in 1:PDE.dim
        fn[dir] .= fn[dir] .+ fp[dir]
    end

    return fn
end

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::Union{AbstractMatrix,Nothing}, numflux::ECChandrashekarNumFlux, PDE::EulerPerfGas)
    return Euler_numflux_Chandrashekar(up, un, PDE)
end

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::AbstractMatrix, numflux::ESChandrashekarDissipNumFlux, PDE::EulerPerfGas)
    f = Euler_numflux_Chandrashekar(up, un, PDE)
    d = Euler_numdissip_ES(un, up, nphys, PDE)

    for dir in 1:PDE.dim
        f[dir] .= f[dir] .+ d[dir]
    end

    return f
end

#####################################################################
# Two-point Fluxes
#####################################################################

abstract type TPFlux end

struct ECSplitTPFlux <: TPFlux end

struct AVSplitTPFlux <: TPFlux end

struct SplitTPFlux <: TPFlux
    alpha::Float64
    beta::Float64
end

struct ECChandrashekarTPFlux <: TPFlux end


# Burgers two-point fluxes

function compute_two_pt_flux(up, un, tpflux::ECSplitTPFlux, PDE::Burgers)
    return ((1/6) .* (un.^2 .+ up .* un .+ up.^2),)
end

function compute_two_pt_flux(up, un, tpflux::AVSplitTPFlux, PDE::Burgers)
    alpha = 0.5
    beta = 0.5
    return (0.25 * alpha .* (un.^2 .+ up.^2) .+ 0.5 * beta .* up .* un,)
end

function compute_two_pt_flux(up, un, tpflux::SplitTPFlux, PDE::Burgers)
    return (0.25 * tpflux.alpha .* (un.^2 .+ up.^2) .+ 0.5 * tpflux.beta .* up .* un,)
end


# Euler two-point fluxes

function compute_two_pt_flux(up, un, tpflux::ECChandrashekarTPFlux, PDE::EulerPerfGas)
    return Euler_numflux_Chandrashekar(up, un, PDE)
end


# Assemble two-point fluxes

function compute_two_pt_flux!(F::Union{Tuple{AbstractArray}, Tuple{AbstractArray, AbstractArray}}, u::AbstractMatrix, uf::AbstractMatrix, tpflux::TPFlux, PDE::GoverningPDE)
    M = size(u,1)
    Npts = size(F[1], 1)

    @inbounds for j in 1:Npts
        if j > M
            up = (@view uf[j - M,:])
        else
            up = (@view u[j,:])
        end
        
        @inbounds for i in 1:(j-1) # Note that we "skip" the diagonal!
            if i > M
                un = (@view uf[i - M,:])
            else
                un = (@view u[i,:])
            end

            f = compute_two_pt_flux(up, un, tpflux, PDE)

            @inbounds for dir in 1:PDE.dim
                @views F[dir][i,j,:] = f[dir]
                @views F[dir][j,i,:] = f[dir]
            end
        end
    end
        
    return F
end


#####################################################################
# Artifical viscosity coefficient
#####################################################################

abstract type ArtViscModel end

struct AVdissip <: ArtViscModel
    addvisc::Union{Float64, Nothing}
end

struct AVEC <: ArtViscModel
    addvisc::Union{Float64, Nothing}
end

struct NoAV <: ArtViscModel
    addvisc::Union{Float64, Nothing}
end

function delta_clip(delta, artviscmodel::AVdissip)
    return -min(0.0, delta)
end

function delta_clip(delta, artviscmodel::AVEC)
    return -delta
end

function delta_clip(delta, artviscmodel::NoAV)
    return 0.0
end

function AV_coeff(delta, den, artviscmodel)
    tol = 1e-14 # tolerance to avoid vanishing denominator

    a = delta_clip(delta, artviscmodel)

    # Clip entropy deficit to avoid error leakage
    if abs(a) < 2.5e-15
        a = 0.0
    end

    if abs(den) < 1e-12 && (abs(delta) / abs(den)) > 1
        println("Careful, visc denominator < 1e-12!")
    end

    if isnothing(artviscmodel.addvisc)
        return a * den / (tol + den^2)
    else
        return a * den / (tol + den^2) + artviscmodel.addvisc
    end
end


#####################################################################
# Residual corrections
#####################################################################

abstract type ResCorrModel end

struct ResCorrEC <: ResCorrModel end

struct ResCorrDissip <: ResCorrModel end

struct NoResCorr <: ResCorrModel end

function delta_clip(delta, rescorr::ResCorrEC)
    return delta
end

function delta_clip(delta, rescorr::ResCorrDissip)
    return min(delta, 0)
end

function delta_clip(delta, rescorr::NoResCorr)
    return 0
end


#####################################################################
# Physics Structure
#####################################################################

# The PhysProp structure contains the PDE, the BCs, the numerical flux and the two-point flux used by the solvers
# The PhysProp structure governs the specific TRANSIENT evolution of the dynamical system.

mutable struct PhysProp
    PDE::GoverningPDE
    source::Union{CLawSource, Nothing}
    BChandler::Dict{Integer, Any}
    numflux::NumFlux
    tpflux::Union{TPFlux,Nothing}
    artvisc::Union{ArtViscModel, Nothing}
    rescorr::Union{ResCorrModel, Nothing}
end


#####################################################################
# Euler helper functions
#####################################################################

# (\rho * e)(u)
function Euler_cvar_intenergy(u::AbstractMatrix, PDE::EulerPerfGas)
    if PDE.dim == 1
        return u[:,end] .- 0.5 .* u[:,2].^2 ./ u[:,1]
    elseif PDE.dim == 2
        return u[:,end] .- 0.5 .* (u[:,2].^2 .+ u[:,3].^2) ./ u[:,1]
    end
end


function Euler_cvar_intenergy(u::AbstractVector, PDE::EulerPerfGas)
    if PDE.dim == 1
        return u[end] - 0.5 * u[2]^2 / u[1]
    elseif PDE.dim == 2
        return u[end] - 0.5 * (u[2]^2 + u[3]^2) / u[1]
    end
end


# (\rho * e)(v)
function Euler_evar_intenergy(v::AbstractMatrix, PDE::EulerPerfGas)
    return ((PDE.gamma-1) ./ (-v[:,end]).^PDE.gamma).^(1/(PDE.gamma-1)) .* exp.(-Euler_evar_entropy(v,PDE) ./ (PDE.gamma-1))
end

function Euler_evar_intenergy(v::AbstractVector, PDE::EulerPerfGas)
    return ((PDE.gamma-1) / (-v[end])^PDE.gamma)^(1/(PDE.gamma-1)) * exp(-Euler_evar_entropy(v,PDE) / (PDE.gamma-1))
end


# (s / cv)(u)
function Euler_cvar_entropy(u::AbstractMatrix, PDE::EulerPerfGas)
    return log.(Euler_pressure(u, PDE) ./ u[:,1].^PDE.gamma)
end

function Euler_cvar_entropy(u::AbstractVector, PDE::EulerPerfGas)
    return log(Euler_pressure(u, PDE) / u[1]^PDE.gamma)
end


# (s / cv)(v)
function Euler_evar_entropy(v::AbstractMatrix, PDE::EulerPerfGas)
    if PDE.dim == 1
        return PDE.gamma .- v[:,1] .+ 0.5 .* v[:,2].^2 ./ v[:,end]
    elseif PDE.dim == 2
        return PDE.gamma .- v[:,1] .+ 0.5 .* (v[:,2].^2 .+ v[:,3].^2) ./ v[:,end]
    end
end

function Euler_evar_entropy(v::AbstractVector, PDE::EulerPerfGas)
    if PDE.dim == 1
        return PDE.gamma - v[1] + 0.5 * v[2]^2 / v[end]
    elseif PDE.dim == 2
        return PDE.gamma - v[1] + 0.5 * (v[2]^2 + v[3]^2) ./ v[end]
    end
end


# v(u)
function Euler_evar(u::AbstractMatrix, PDE::EulerPerfGas)
    rhoe = Euler_cvar_intenergy(u,PDE)
    s = Euler_cvar_entropy(u, PDE)

    v = Array{eltype(u)}(undef, size(u)...)
    @views v[:,1] .= (-s .+ PDE.gamma .+ 1) .- u[:,end] ./ rhoe
    @views v[:,2:1+PDE.dim] .= u[:,2:1+PDE.dim] ./ rhoe
    @views v[:,end] .= -u[:,1] ./ rhoe

    return v
end

function Euler_evar(u::AbstractVector, PDE::EulerPerfGas)
    rhoe = Euler_cvar_intenergy(u,PDE)
    s = Euler_cvar_entropy(u, PDE)

    v = Vector{eltype(u)}(undef, PDE.dim+2)
    v[1] = (-s + PDE.gamma + 1) - u[end] / rhoe
    v[2:1+PDE.dim] .= u[2:1+PDE.dim] ./ rhoe
    v[end] = -u[1] / rhoe

    return v
end


# u(v)
function Euler_cvar(v::AbstractMatrix, PDE::EulerPerfGas)
    rhoe = Euler_evar_intenergy(v, PDE)

    u = Array{eltype(v)}(undef, size(v)...)
    @views u[:,1] .= -rhoe .* v[:,end]
    @views u[:,2:1+PDE.dim] .= v[:,2:1+PDE.dim] .* rhoe
    @views u[:,end] .= rhoe .* (1 .- 0.5 .* sum(v[:,2:1+PDE.dim].^2, dims=2) ./  v[:,end])

    return u
end

function Euler_cvar(v::AbstractVector, PDE::EulerPerfGas)
    rhoe = Euler_evar_intenergy(v, PDE)

    u = Vector{eltype(v)}(undef, PDE.dim+2)
    u[1] = -rhoe * v[end]
    @views u[2:1+PDE.dim] .= v[2:1+PDE.dim] .* rhoe
    u[end] = rhoe * (1 - 0.5 * sum(v[2:1+PDE.dim].^2) /  v[end])

    return u
end


# p(u)
function Euler_pressure(u::AbstractMatrix, PDE::EulerPerfGas)
    return (PDE.gamma - 1) .* Euler_cvar_intenergy(u, PDE)
end

function Euler_pressure(u::AbstractVector, PDE::EulerPerfGas)
    return (PDE.gamma - 1) * Euler_cvar_intenergy(u, PDE)
end


# (\partial u / \partial v)(u)
function Euler_cvar_Hessian(u::AbstractVector, PDE::EulerPerfGas) # copied from (Chan, 2025)
    if PDE.dim == 1
        K = Matrix{eltype(u)}(undef,3,3)

        p = Euler_pressure(u, PDE)
        a2 = PDE.gamma * Euler_pressure(u, PDE) / u[1]

        K[1,:] = u
        K[2,2] = u[2]^2 / u[1] + p
        K[2,3] = u[2] / u[1] * (u[end] + p)
        K[3,3] = u[1] * (a2 / (PDE.gamma-1) + 0.5 * u[2]^2 / u[1])^2 - a2 * p / (PDE.gamma - 1)

        copyto!(K, Symmetric(K, :U)) # symmetrize

        return K
    else
        error("Euler Hessian only implemented in 1D.")
    end
end


# logmean
function logmean(up::Real, un::Real)
    tol = 1e-2 # tolerance for logmean computation

    if abs(up - un) < tol # Roe simplification for up -> un
        r = up/un
        a = ((r - 1) / (r + 1))^2
        return 0.5 * (up + un) / (1 + a/3 + a^2/5 + a^3/7)
    else
        return (up - un) / (log(up) - log(un))
    end
end


# f_dir
function Euler_physflux(u::AbstractMatrix, PDE::EulerPerfGas)
    p = Euler_pressure(u, PDE)

    f1 = Array{eltype(u)}(undef, size(u)...)
    @views f1[:,1] .= u[:,2]
    @views f1[:,2:PDE.dim+1] .=  u[:,2:PDE.dim+1] .* u[:,2] ./  u[:,1]
    @views f1[:,2] .= f1[:,2] .+ p
    @views f1[:,end] .= (p .+  u[:,end]) .* (u[:,2]./ u[:,1])

    if PDE.dim == 1
        return (f1,)
    elseif PDE.dim == 2
        f2 = Array{eltype(u)}(undef, size(u)...)
        @views f2[:,1] .= u[3,:]
        @views f2[:,2] .= f1[:,3]
        @views f2[:,3] .= u[:,3].^2 ./  u[:,1] .+ p
        @views f2[:,end] .= (p .+  u[:,end]) .* (u[:,3]./ u[:,1])

        return (f1, f2)
    end
end

function Euler_physflux(u::AbstractVector, PDE::EulerPerfGas)
    p = Euler_pressure(u, PDE)

    f1 = Vector{eltype(u)}(undef, PDE.dim+2)
    f1[1] = u[2]
    @views f1[2:PDE.dim+1] .=  u[2:PDE.dim+1] .* u[2] ./  u[1]
    f1[2] = f1[2] + p
    f1[end] = (p + u[end]) * u[2] / u[1]

    if PDE.dim == 1
        return (f1,)
    elseif PDE.dim == 2
        f2 = Vector{eltype(u)}(undef, PDE.dim+2)
        f2[1] = u[3]
        f2[2] = f1[3]
        f2[3] = u[3]^2 /  u[1] + p
        f2[end] = (p + u[end]) * (u[3]/ u[1])

        return (f1, f2)
    end

end


# Numerical fluxes, ONLY FOR 1D right now
function Euler_numflux_Chandrashekar(un::AbstractMatrix, up::AbstractMatrix, PDE::EulerPerfGas) # copied from (Chan 2018)
    if PDE.dim == 1
        f1 = Array{eltype(up)}(undef, size(up)...)

        @views velp = up[:,2:PDE.dim+1] ./ up[:,1]
        @views veln = un[:,2:PDE.dim+1] ./ un[:,1]
        velavg = 0.5 .* (velp .+ veln)

        @views betap = 0.5 .* up[:,1] ./ Euler_pressure(up, PDE)
        @views betan = 0.5 .* un[:,1] ./ Euler_pressure(un, PDE)

        @views @. f1[:,1] = logmean(up[:,1], un[:,1]) * velavg[:,1]
        @views @. f1[:,2] = 0.5 * (up[:,1] + un[:,1]) / (betap + betan) + velavg * f1[:,1]
        @views @. f1[:,3] = f1[:,1] * (0.5 / (PDE.gamma-1) / logmean(betan,betap) - 0.25 * (velp^2 + veln^2)) + velavg * f1[:,2]

        return (f1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end

function Euler_numflux_Chandrashekar(un::AbstractVector, up::AbstractVector, PDE::EulerPerfGas) # copied from (Chan 2018)
    if PDE.dim == 1
        f1 = Vector{eltype(up)}(undef, 3)

        velp = up[2] / up[1]
        veln = un[2] / un[1]
        velavg = 0.5 * (velp + veln)

        betap = 0.5 * up[1] / Euler_pressure(up, PDE)
        betan = 0.5 * un[1] / Euler_pressure(un, PDE)

        f1[1] = logmean(up[1], un[1]) * velavg[1]
        f1[2] = 0.5 * (up[1] + un[1]) / (betap + betan) + velavg * f1[1]
        f1[3] = f1[1] * (0.5 / (PDE.gamma-1) / logmean(betan,betap) - 0.25 * (velp^2 + veln^2)) + velavg * f1[2]

        return (f1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end


function Euler_numdissip_ES(un::AbstractMatrix, up::AbstractMatrix, nphys::AbstractMatrix, PDE::EulerPerfGas) # copied from (Gassner, Kopriva, Winters, Hindenlang, 2018)
    if PDE.dim == 1
        d1 = Array{eltype(up)}(undef, size(up)...)

        wjump = (Euler_evar(up, PDE) - Euler_evar(un, PDE)) .* nphys # WOULD BE MUCH MORE EFFICIENT IF WE COULD USE THE ENTROPY VARS AS INPUT

        # We iterate to find the diffusion term (too expensive to do in omne shot)
        R = [1.0 1.0 1.0; 0 0 0 ; 0 0 0] # prealloc for R

        for index in axes(d1,1)
            velp = up[index,2] ./ up[index,1]
            veln = un[index,2] ./ un[index,1]
            velavg = 0.5 * (velp + veln)
            vel2avg = 2 * velavg^2 - 0.5 * (veln^2 + veln^2)

            pn = Euler_pressure(up[index,:], PDE)
            pp = Euler_pressure(un[index,:], PDE)

            rholn = logmean(up[index,1], un[index,1])
            betaln = logmean(0.5 * up[index,1] / pn, 0.5 * un[index,1] / pp)

            abar = sqrt(0.5 * PDE.gamma * (pn + pp)/rholn)
            hbar = PDE.gamma/(2.0*betaln*(PDE.gamma-1.0)) + 0.5 * vel2avg

            R[2,1] = velavg - abar ; R[2,2] = velavg ; R[2,3] = velavg + abar
            R[3,1] = hbar - velavg*abar ; R[3,2] = 0.5*vel2avg ; R[3,3] = hbar + velavg * abar

            Lambda = abs.([velavg - abar, velavg, velavg + abar])
            T = [rholn/2.0/PDE.gamma, rholn * (PDE.gamma-1.0)/PDE.gamma, rholn/2.0/PDE.gamma]

            d1[index,:] = -0.5 * R * Diagonal(Lambda) * Diagonal(T) * R' * wjump[index,:]
        end

        return (d1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end


function Euler_numdissip_ES(un::AbstractVector, up::AbstractVector, nphys::AbstractVector, PDE::EulerPerfGas) # copied from (Chan 2018)
    if PDE.dim == 1

        wjump = (Euler_evar(up, PDE) - Euler_evar(un, PDE)) .* nphys[1] # WOULD BE MUCH MORE EFFICIENT IF WE COULD USE THE ENTROPY VARS AS INPUT

        velp = up[2] ./ up[1]
        veln = un[2] ./ un[1]
        velavg = 0.5 * (velp + veln)
        vel2avg = 2 * velavg^2 - 0.5 * (veln^2 + veln^2)

        pn = Euler_pressure(up, PDE)
        pp = Euler_pressure(un, PDE)

        rholn = logmean(up[1], un[1])
        betaln = logmean(0.5 * up[1] / pn, 0.5 * un[1] / pp)

        abar = sqrt(0.5 * PDE.gamma * (pn + pp)/rholn)
        hbar = PDE.gamma/(2.0*betaln*(PDE.gamma-1.0)) + 0.5 * vel2avg

        R = [1.0 1.0 1.0;
             velavg - abar velavg velavg + abar ;
             hbar - velavg*abar 0.5*vel2avg hbar + velavg * abar]

        Lambda = abs.([velavg - abar, velavg, velavg + abar])
        T = [rholn/2.0/PDE.gamma, rholn * (PDE.gamma-1.0)/PDE.gamma, rholn/2.0/PDE.gamma]

        d1 = -0.5 * R * Diagonal(Lambda) * Diagonal(T) * R' * wjump

        return (d1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end
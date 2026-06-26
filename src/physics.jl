#####################################################################
# Everything concerning the problem's physics
#####################################################################

# THE SOLUTION IS ALWAYS A MATRIX OF SIZE (Nstate, NDOF)

#####################################################################
# Physical Fluxes
#####################################################################

function compute_physflux(u::AbstractMatrix, param::parameters)
    if param.pdetype == "LinAdv"
        return (param.a .* u,)

    elseif param.pdetype == "Burgers"
        return (0.5 .* u.^2,)

    elseif param.pdetype == "EulerPerfGas"
        return Euler_physflux(u, param)
    else
        error("Unknown PDE type!")
    end
end

#####################################################################
# Numerical Fluxes
#####################################################################

function compute_numflux(un::AbstractMatrix, up::AbstractMatrix, nphys::Union{AbstractMatrix,Nothing}, param::parameters)
    if param.pdetype == "LinAdv"
        if param.numfluxtype == "central"
            return (0.5 .* param.a .* (up .+ un),)
        elseif param.numfluxtype == "upwind"
            f = copy(up)
            if param.a * nphys[1] >= 0
                index = nphys .== nphys[1]
            else
                index = nphys .!= nphys[1]
            end
            f[index] = un[index]  
            return (param.a .* f,)
        else
            error("Undefined numerical flux!")
        end

    elseif param.pdetype == "Burgers"
        if param.numfluxtype == "central"
            return (0.25 .* (up.^2 .+ un.^2),)
        elseif param.numfluxtype == "LF"
            return (0.25 .* ((un.^2 .+ up.^2) .- max.(abs.(up), abs.(un)) .* (up .- un) .* nphys),)
        elseif param.numfluxtype == "upwind"
            f = similar(up)
            vel = 0.5 .* (un.^2 .- up.^2) ./ (un .- up) # local velocity
            slicing = vel .* nphys .>= 0
            notslicing = .!slicing
            f[slicing] .= 0.5 .* un[slicing].^2
            f[notslicing] .= 0.5 .* up[notslicing].^2
            return (f,)
        elseif param.numfluxtype == "EC_split"
            return ((1/6) .* (un.^2 .+ up .* un .+ up.^2),)
        else
            error("Undefined numerical flux!")
        end

    elseif param.pdetype == "EulerPerfGas"
        if param.numfluxtype == "central"
            fp = Euler_physflux(up, param)
            fn = Euler_physflux(up, param)
            for dir in 1:param.dim
                fn[dir] .= fn[dir] .+ fp[dir]
            end
            return fn

        elseif param.numfluxtype == "EC_Chandrashekar"
            return Euler_numflux_Chandrashekar(up, un, param)

        elseif param.numfluxtype == "ES_Chandrashekar_dissip"
            f = Euler_numflux_Chandrashekar(up, un, param)
            d = Euler_numdissip_ES(un, up, nphys, param)
            for dir in 1:param.dim
                f[dir] .= f[dir] .+ d[dir]
            end
            return f

        else
            error("Undefined numerical flux!")
        end
    end
end

#####################################################################
# Two-point Fluxes
#####################################################################

function compute_two_pt_flux!(F::Union{Tuple{AbstractArray}, Tuple{AbstractArray, AbstractArray}}, u::AbstractMatrix, uf::AbstractMatrix, param::parameters)
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

            if param.pdetype == "Burgers"
                if param.twoptfluxtype == "EC_split"
                    f = ((1/6) .* (un.^2 .+ up .* un .+ up.^2),)
                elseif param.twoptfluxtype == "AV_split"
                    f = ((1/8) .* (un.^2 .+ 2.0 .* up .* un .+ up.^2),)
                elseif param.twoptfluxtype == "OQ_split"
                    p = M-1
                    f = ((1/(4.0 * (2*p+1))) .* ((p+1).*(un.^2 .+ up.^2) .+ (2.0 * p) .* (up .* un)),)
                else
                    error("Undefined two-point flux!")
                end
            elseif param.pdetype == "EulerPerfGas"
                if param.twoptfluxtype == "EC_Chandrashekar"
                    f = Euler_numflux_Chandrashekar(up, un, param)
                else
                    error("Undefined two-point flux!")
                end
            else
                error("Unknown PDE type!")
            end

            @inbounds for dir in 1:param.dim
                @views F[dir][i,j,:] = f[dir]
                @views F[dir][j,i,:] = f[dir]
            end
        end
    end
        
    return F
end

#####################################################################
# Entropy mappings
#####################################################################

function compute_evar(u::AbstractMatrix, param::parameters)
    if param.pdetype == "Burgers"
        return u
    elseif param.pdetype == "EulerPerfGas"
        return Euler_evar(u, param)
    end
end

function compute_cvar(v::AbstractMatrix, param::parameters)
    if param.pdetype == "Burgers"
        return v
    elseif param.pdetype == "EulerPerfGas"
        return Euler_cvar(v, param)
    end
end

#####################################################################
# Compute Entropy
#####################################################################

function compute_local_entropy(u::AbstractMatrix, param::parameters)
    if param.pdetype == "LinAdv"
        return 0.5 .* u.^2
    elseif param.pdetype == "Burgers"
        return 0.5 .* u.^2
    elseif param.pdetype == "EulerPerfGas"
        return -u[:,1] .* Euler_cvar_entropy(u, param)
    end
end

#####################################################################
# Compute Entropy Potentials
#####################################################################

function compute_cvar_potential(u::AbstractArray, param::parameters)
    if param.pdetype == "Burgers"
        return ((1/6) .* u.^3,)
    elseif param.pdetype == "EulerPerfGas"
        return Tuple((param.gamma-1) .* u[:,istate] for istate in 2:1+param.dim) # CAREFUL, MISTAKE IN (Chan, 2025)
    end
end

#####################################################################
# Compute Entropy Hessian
#####################################################################

function compute_cvar_Hessian(u::AbstractVector, param::parameters) # Operates on one node at the time.
    if param.pdetype == "Burgers"
        K = Matrix{Float64}(undef, 1, 1)
        K[1,1] = 1.0
        return K
    elseif param.pdetype == "EulerPerfGas"
        return Euler_cvar_Hessian(u, param)
    end
end

#####################################################################
# Artifical viscosity coefficient
#####################################################################

function AV_coeff(delta, den, param)
    tol = 1e-14 # tolerance to avoid vanishing denominator

    if param.AVcoeff == "AVdissip"
        a = -min(0.0, delta)
    elseif param.AVcoeff == "AVEC"
        a = -delta
    elseif param.AVcoeff == "NoAV"
        a = 0.0
    end

    # Clip entropy deficit to avoid error leakage
    if abs(a) < 2.5e-15
        a = 0.0
    end

    if abs(den) < 1e-12 && (abs(delta) / abs(den)) > 1
        println("Careful, visc denominator < 1e-12!")
    end

    if isnothing(param.addviscosity)
        return a * den / (tol + den^2)
    else
        return a * den / (tol + den^2) + param.addviscosity
    end

end

#####################################################################
# Euler helper functions
#####################################################################

# (\rho * e)(u)
function Euler_cvar_intenergy(u::AbstractMatrix, param::parameters)
    if param.dim == 1
        return u[:,end] .- 0.5 .* u[:,2].^2 ./ u[:,1]
    elseif param.dim == 2
        return u[:,end] .- 0.5 .* (u[:,2].^2 .+ u[3,:].^2) ./ u[1,:]
    end
end


function Euler_cvar_intenergy(u::AbstractVector, param::parameters)
    if param.dim == 1
        return u[end] - 0.5 * u[2]^2 / u[1]
    elseif param.dim == 2
        return u[end] - 0.5 * (u[2]^2 + u[3]^2) / u[1]
    end
end


# (\rho * e)(v)
function Euler_evar_intenergy(v::AbstractMatrix, param::parameters)
    return ((param.gamma-1) ./ (-v[:,end]).^param.gamma).^(1/(param.gamma-1)) .* exp.(-Euler_evar_entropy(v,param) ./ (param.gamma-1))
end

function Euler_evar_intenergy(v::AbstractVector, param::parameters)
    return ((param.gamma-1) / (-v[end])^param.gamma)^(1/(param.gamma-1)) * exp(-Euler_evar_entropy(v,param) / (param.gamma-1))
end


# (s / cv)(u)
function Euler_cvar_entropy(u::AbstractMatrix, param::parameters)
    return log.(Euler_pressure(u, param) ./ u[:,1].^param.gamma)
end

function Euler_cvar_entropy(u::AbstractVector, param::parameters)
    return log(Euler_pressure(u, param) / u[1]^param.gamma)
end


# (s / cv)(v)
function Euler_evar_entropy(v::AbstractMatrix, param::parameters)
    if param.dim == 1
        return param.gamma .- v[:,1] .+ 0.5 .* v[:,2].^2 ./ v[:,end]
    elseif param.dim == 2
        return param.gamma .- v[:,1] .+ 0.5 .* (v[:,2].^2 .+ v[:,3].^2) ./ v[:,end]
    end
end

function Euler_evar_entropy(v::AbstractVector, param::parameters)
    if param.dim == 1
        return param.gamma - v[1] + 0.5 * v[2]^2 / v[end]
    elseif param.dim == 2
        return param.gamma - v[1] + 0.5 * (v[2]^2 + v[3]^2) ./ v[end]
    end
end


# v(u)
function Euler_evar(u::AbstractMatrix, param::parameters)
    rhoe = Euler_cvar_intenergy(u,param)
    s = Euler_cvar_entropy(u, param)

    v = Array{eltype(u)}(undef, size(u)...)
    @views v[:,1] .= (-s .+ param.gamma .+ 1) .- u[:,end] ./ rhoe
    @views v[:,2:1+param.dim] .= u[:,2:1+param.dim] ./ rhoe
    @views v[:,end] .= -u[:,1] ./ rhoe

    return v
end

function Euler_evar(u::AbstractVector, param::parameters)
    rhoe = Euler_cvar_intenergy(u,param)
    s = Euler_cvar_entropy(u, param)

    v = Vector{eltype(u)}(undef, param.dim+2)
    v[1] = (-s + param.gamma + 1) - u[end] / rhoe
    v[2:1+param.dim] .= u[2:1+param.dim] ./ rhoe
    v[end] = -u[1] / rhoe

    return v
end


# u(v)
function Euler_cvar(v::AbstractMatrix, param::parameters)
    rhoe = Euler_evar_intenergy(v, param)

    u = Array{eltype(v)}(undef, size(v)...)
    @views u[:,1] .= -rhoe .* v[:,end]
    @views u[:,2:1+param.dim] .= v[:,2:1+param.dim] .* rhoe
    @views u[:,end] .= rhoe .* (1 .- 0.5 .* sum(v[:,2:1+param.dim].^2, dims=2) ./  v[:,end])

    return u
end

function Euler_cvar(v::AbstractVector, param::parameters)
    rhoe = Euler_evar_intenergy(v, param)

    u = Vector{eltype(v)}(undef, param.dim+2)
    u[1] = -rhoe * v[end]
    @views u[2:1+param.dim] .= v[2:1+param.dim] .* rhoe
    u[end] = rhoe * (1 - 0.5 * sum(v[2:1+param.dim].^2) /  v[end])

    return u
end


# p(u)
function Euler_pressure(u::AbstractMatrix, param::parameters)
    return (param.gamma - 1) .* Euler_cvar_intenergy(u, param)
end

function Euler_pressure(u::AbstractVector, param::parameters)
    return (param.gamma - 1) * Euler_cvar_intenergy(u, param)
end


# (\partial u / \partial v)(u)
function Euler_cvar_Hessian(u::AbstractVector, param::parameters) # copied from (Chan, 2025)
    if param.dim == 1
        K = Matrix{eltype(u)}(undef,3,3)

        p = Euler_pressure(u, param)
        a2 = param.gamma * Euler_pressure(u, param) / u[1]

        K[1,:] = u
        K[2,2] = u[2]^2 / u[1] + p
        K[2,3] = u[2] / u[1] * (u[end] + p)
        K[3,3] = u[1] * (a2 / (param.gamma-1) + 0.5 * u[2]^2 / u[1])^2 - a2 * p / (param.gamma - 1)

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
function Euler_physflux(u::AbstractMatrix, param::parameters)
    p = Euler_pressure(u, param)

    f1 = Array{eltype(u)}(undef, size(u)...)
    @views f1[:,1] .= u[:,2]
    @views f1[:,2:param.dim+1] .=  u[:,2:param.dim+1] .* u[:,2] ./  u[:,1]
    @views f1[:,2] .= f1[:,2] .+ p
    @views f1[:,end] .= (p .+  u[:,end]) .* (u[:,2]./ u[:,1])

    if param.dim == 1
        return (f1,)
    elseif param.dim == 2
        f2 = Array{eltype(u)}(undef, size(u)...)
        @views f2[:,1] .= u[3,:]
        @views f2[:,2] .= f1[:,3]
        @views f2[:,3] .= u[:,3].^2 ./  u[:,1] .+ p
        @views f2[:,end] .= (p .+  u[:,end]) .* (u[:,3]./ u[:,1])

        return (f1, f2)
    end
end

function Euler_physflux(u::AbstractVector, param::parameters)
    p = Euler_pressure(u, param)

    f1 = Vector{eltype(u)}(undef, param.dim+2)
    f1[1] = u[2]
    @views f1[2:param.dim+1] .=  u[2:param.dim+1] .* u[2] ./  u[1]
    f1[2] = f1[2] + p
    f1[end] = (p + u[end]) * u[2] / u[1]

    if param.dim == 1
        return (f1,)
    elseif param.dim == 2
        f2 = Vector{eltype(u)}(undef, param.dim+2)
        f2[1] = u[3]
        f2[2] = f1[3]
        f2[3] = u[3]^2 /  u[1] + p
        f2[end] = (p + u[end]) * (u[3]/ u[1])

        return (f1, f2)
    end

end


# Numerical fluxes, ONLY FOR 1D right now
function Euler_numflux_Chandrashekar(un::AbstractMatrix, up::AbstractMatrix, param::parameters) # copied from (Chan 2018)
    if param.dim == 1
        f1 = Array{eltype(up)}(undef, size(up)...)

        @views velp = up[:,2:param.dim+1] ./ up[:,1]
        @views veln = un[:,2:param.dim+1] ./ un[:,1]
        velavg = 0.5 .* (velp .+ veln)

        @views betap = 0.5 .* up[:,1] ./ Euler_pressure(up, param)
        @views betan = 0.5 .* un[:,1] ./ Euler_pressure(un, param)

        @views @. f1[:,1] = logmean(up[:,1], un[:,1]) * velavg[:,1]
        @views @. f1[:,2] = 0.5 * (up[:,1] + un[:,1]) / (betap + betan) + velavg * f1[:,1]
        @views @. f1[:,3] = f1[:,1] * (0.5 / (param.gamma-1) / logmean(betan,betap) - 0.25 * (velp^2 + veln^2)) + velavg * f1[:,2]

        return (f1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end

function Euler_numflux_Chandrashekar(un::AbstractVector, up::AbstractVector, param::parameters)
    if param.dim == 1
        f1 = Vector{eltype(up)}(undef, 3)

        velp = up[2] / up[1]
        veln = un[2] / un[1]
        velavg = 0.5 * (velp + veln)

        betap = 0.5 * up[1] / Euler_pressure(up, param)
        betan = 0.5 * un[1] / Euler_pressure(un, param)

        f1[1] = logmean(up[1], un[1]) * velavg[1]
        f1[2] = 0.5 * (up[1] + un[1]) / (betap + betan) + velavg * f1[1]
        f1[3] = f1[1] * (0.5 / (param.gamma-1) / logmean(betan,betap) - 0.25 * (velp^2 + veln^2)) + velavg * f1[2]

        return (f1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end


function Euler_numdissip_ES(un::AbstractMatrix, up::AbstractArray, nphys::AbstractMatrix, param::parameters) # copied from (Gassner, Kopriva, Winters, Hindenlang, 2018)
    if param.dim == 1
        d1 = Array{eltype(up)}(undef, size(up)...)

        wjump = (Euler_evar(up, param) - Euler_evar(un, param)) .* nphys # WOULD BE MUCH MORE EFFICIENT IF WE COULD USE THE ENTROPY VARS AS INPUT

        # We iterate to find the diffusion term (too expensive to do in omne shot)
        R = [1.0 1.0 1.0; 0 0 0 ; 0 0 0] # prealloc for R

        for index in axes(d1,1)
            velp = up[index,2] ./ up[index,1]
            veln = un[index,2] ./ un[index,1]
            velavg = 0.5 * (velp + veln)
            vel2avg = 2 * velavg^2 - 0.5 * (veln^2 + veln^2)

            pn = Euler_pressure(up[index,:], param)
            pp = Euler_pressure(un[index,:], param)

            rholn = logmean(up[index,1], un[index,1])
            betaln = logmean(0.5 * up[index,1] / pn, 0.5 * un[index,1] / pp)

            abar = sqrt(0.5 * param.gamma * (pn + pp)/rholn)
            hbar = param.gamma/(2.0*betaln*(param.gamma-1.0)) + 0.5 * vel2avg

            R[2,1] = velavg - abar ; R[2,2] = velavg ; R[2,3] = velavg + abar
            R[3,1] = hbar - velavg*abar ; R[3,2] = 0.5*vel2avg ; R[3,3] = hbar + velavg * abar

            Lambda = abs.([velavg - abar, velavg, velavg + abar])
            T = [rholn/2.0/param.gamma, rholn * (param.gamma-1.0)/param.gamma, rholn/2.0/param.gamma]

            d1[index,:] = -0.5 * R * Diagonal(Lambda) * Diagonal(T) * R' * wjump[index,:]
        end

        return (d1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end

function Euler_numdissip_ES(un::AbstractVector, up::AbstractVector, nphys::AbstractVector, param::parameters) # copied from (Chan 2018)
    if param.dim == 1

        wjump = (Euler_evar(up, param) - Euler_evar(un, param)) .* nphys[1] # WOULD BE MUCH MORE EFFICIENT IF WE COULD USE THE ENTROPY VARS AS INPUT

        velp = up[2] ./ up[1]
        veln = un[2] ./ un[1]
        velavg = 0.5 * (velp + veln)
        vel2avg = 2 * velavg^2 - 0.5 * (veln^2 + veln^2)

        pn = Euler_pressure(up, param)
        pp = Euler_pressure(un, param)

        rholn = logmean(up[1], un[1])
        betaln = logmean(0.5 * up[1] / pn, 0.5 * un[1] / pp)

        abar = sqrt(0.5 * param.gamma * (pn + pp)/rholn)
        hbar = param.gamma/(2.0*betaln*(param.gamma-1.0)) + 0.5 * vel2avg

        R = [1.0 1.0 1.0;
             velavg - abar velavg velavg + abar ;
             hbar - velavg*abar 0.5*vel2avg hbar + velavg * abar]

        Lambda = abs.([velavg - abar, velavg, velavg + abar])
        T = [rholn/2.0/param.gamma, rholn * (param.gamma-1.0)/param.gamma, rholn/2.0/param.gamma]

        d1 = -0.5 * R * Diagonal(Lambda) * Diagonal(T) * R' * wjump

        return (d1,)
    else
        error("Chandrashekar flux only implemented in 1D.")
    end
end
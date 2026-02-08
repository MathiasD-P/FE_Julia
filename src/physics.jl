#####################################################################
# Everything concerning the problem's physics
#####################################################################

# THE SOLUTION IS ALWAYS A MATRIX OF SIZE (Nstate, NDOF)

#####################################################################
# Physical Fluxes
#####################################################################

function compute_physflux(u::AbstractArray, param::parameters)
    if param.pdetype == "LinAdv"
        return [param.a .* u]

    elseif param.pdetype == "Burgers"
        return [0.5 .* u.^2]

    elseif param.pdetype == "EulerPerfGas"
        return Euler_physflux(u, param)
    else
        error("Unknown PDE type!")
    end
end

#####################################################################
# Numerical Fluxes
#####################################################################

function compute_numflux(un::AbstractArray, up::AbstractArray, nphys::Union{AbstractArray,Nothing}, param::parameters)
    if param.pdetype == "LinAdv"
        if param.numfluxtype == "central"
            return [0.5 .* param.a .* (up .+ un)]
        elseif param.numfluxtype == "upwind"
            f = copy(up)
            if param.a * nphys[1] >= 0
                index = nphys .== nphys[1]
            else
                index = nphys .!= nphys[1]
            end
            f[index] = un[index]  
            return [param.a .* f]
        else
            error("Undefined numerical flux!")
        end

    elseif param.pdetype == "Burgers"
        if param.numfluxtype == "central"
            return [0.25 .* (up.^2 .+ un.^2)]
        elseif param.numfluxtype == "LF"
            return [0.25 .* ((un.^2 .+ up.^2) .- max.(abs.(up), abs.(un)) .* (up .- un) .* dg.nphys)]
        elseif param.numfluxtype == "EC_split"
            return [(1/6) .* (un.^2 .+ up .* un .+ up.^2)]
        else
            error("Undefined numerical flux!")
        end

    elseif param.pdetype == "EulerPerfGas"
        if param.numfluxtype == "central"
            return 0.5 .* (Euler_physflux(up, param) .+ Euler_physflux(un, param))
        elseif param.numfluxtype == "EC_Chandrashekar"
            return Euler_numflux_Chandrashekar(up, un, param)
        else
            error("Undefined numerical flux!")
        end
    end
end

#####################################################################
# Two-point Fluxes
#####################################################################

function two_pt_flux(u::AbstractArray, uf::AbstractArray, param::parameters)
    M, Nstates = size(u)
    Npts = M + size(uf, 1)

    F = [Array{Float64}(undef, Npts,Npts,Nstates)]
    if param.dim == 2
        F = [Array{Float64}(undef, Npts,Npts,Nstates), Array{Float64}(undef, Npts,Npts,Nstates)]
    end

    if param.dim == 1
        @inbounds for j in 1:Npts
            @inbounds for i in 1:j
                if j > M
                    up = (@view uf[j - M,:])'
                else
                    up = (@view u[j,:])'
                end
                if i > M
                    un = (@view uf[i - M,:])'
                else
                    un = (@view u[i,:])'
                end

                if param.pdetype == "Burgers"
                    if param.twoptfluxtype == "EC_split"
                        f = [(1/6) .* (un.^2 .+ up .* un .+ up.^2)]
                    else
                        error("Undefined numerical flux!")
                    end
                elseif param.pdetype == "EulerPerfGas"
                    if param.twoptfluxtype == "EC_Chandrashekar"
                        f = Euler_numflux_Chandrashekar(up, un, param)
                    else
                        error("Undefined numerical flux!")
                    end
                else
                    error("Unknown PDE type!")
                end

                if param.dim == 1
                    F[1][i,j,:] .= f[1][:]
                    F[1][j,i,:] .= f[1][:]
                elseif param.dim == 2
                    F[1][i,j,:] .= f[1][:]
                    F[1][j,i,:] .= f[1][:]
                    F[2][i,j,:] .= f[2][:]
                    F[2][j,i,:] .= f[2][:]
                end
            end
        end
        
        return F
    end
end

#####################################################################
# Entropy mappings
#####################################################################

function compute_evar(u::AbstractArray, param::parameters)
    if param.pdetype == "Burgers"
        return u
    elseif param.pdetype == "EulerPerfGas"
        rhoe = Euler_cvar_intenergy(u,param)
        s = Euler_cvar_entropy(u, param)

        v = zeros(size(u))
        @views v[:,1] .= (-s .+ param.gamma .+ 1) .- u[:,end] ./ rhoe
        @views v[:,2:1+param.dim] .= u[:,2:1+param.dim] ./ rhoe
        @views v[:,end] .= -u[:,1] ./ rhoe

        return v
    end
end

function compute_cvar(v::AbstractArray, param::parameters)
    if param.pdetype == "Burgers"
        return v
    elseif param.pdetype == "EulerPerfGas"
        rhoe = Euler_evar_intenergy(v, param)

        u = zeros(size(v))
        @views u[:,1] .= -rhoe .* v[:,end]
        @views u[:,2:1+param.dim] .= v[:,2:1+param.dim] .* rhoe
        @views u[:,end] .= rhoe .* (1 .- 0.5 .* sum(v[:,2:1+param.dim].^2, dims=2) ./  v[:,end])

        return u
    end
end

#####################################################################
# Compute Entropy
#####################################################################

function compute_local_entropy(u::AbstractArray, param::parameters)
    if param.pdetype == "LinAdv"
        return 0.5 .* u.^2
    elseif param.pdetype == "Burgers"
        return 0.5 .* u.^2
    elseif param.pdetype == "EulerPerfGas"
        return -u[:,1] .* Euler_cvar_entropy(u, param)
    end
end

#####################################################################
# Euler helper functions
#####################################################################

# (\rho * e)(u)
function Euler_cvar_intenergy(u::AbstractArray, param::parameters)
    if param.dim == 1
        return u[:,end] .- 0.5 .* u[:,2].^2 ./ u[:,1]
    elseif param.dim == 2
        return u[:,end] .- 0.5 .* (u[:,2].^2 .+ u[3,:].^2) ./ u[1,:]
    end
end

# (\rho * e)(v)
function Euler_evar_intenergy(v::AbstractArray, param::parameters)
    return ((param.gamma-1) ./ (-v[:,end]).^param.gamma).^(1/(param.gamma-1)) .* exp.(-Euler_evar_entropy(v,param) ./ (param.gamma-1))
end

# (s / cv)(u)
function Euler_cvar_entropy(u::AbstractArray, param::parameters)
    return log.(Euler_pressure(u, param) ./ u[:,1].^param.gamma)
end

# (s / cv)(v)
function Euler_evar_entropy(v::AbstractArray, param::parameters)
    if param.dim == 1
        return param.gamma .- v[:,1] .+ 0.5 .* v[:,2].^2 ./ v[:,end]
    elseif param.dim == 2
        return param.gamma .- v[:,1] .+ 0.5 .* (v[:,2].^2 .+ v[:,3].^2) ./ v[:,end]
    end
end

# p(u)
function Euler_pressure(u::AbstractArray, param::parameters)
    return (param.gamma - 1) .* Euler_cvar_intenergy(u, param)
end

# logmean
function logmean(up::Float64, un::Float64)
    epsilon = 1e-2 # tolerance for logmean computation

    if abs(up - un) < epsilon # Roe simplification for up -> un
        r = up/un
        a = ((r - 1) / (r + 1))^2
        return 0.5 * (up + un) / (1 + a/3 + a^2/5 + a^3/7)
    else
        return (up - un) / (log(up) - log(un))
    end
end

# f_dir
function Euler_physflux(u::AbstractArray, param::parameters)
    p = Euler_pressure(u, param)

    f1 = zeros(size(u))
    @views f1[:,1] .= u[:,2]
    @views f1[:,2:param.dim+1] .=  u[:,2:param.dim+1] .* u[:,2] ./  u[:,1]
    @views f1[:,2] .= f1[:,2] .+ p
    @views f1[:,end] .= (p .+  u[:,end]) .* (u[:,2]./ u[:,1])

    if param.dim == 1
        return [f1]
    elseif param.dim == 2
        f2 = zeros(size(u))
        @views f2[:,1] .= u[3,:]
        @views f2[:,2] .= f1[:,3]
        @views f2[:,3] .= u[:,3].^2 ./  u[:,1] .+ p
        @views f2[:,end] .= (p .+  u[:,end]) .* (u[:,3]./ u[:,1])

        return [f1, f2]
    end

end

# ONLY FOR 1D RIGHT NOW
function Euler_numflux_Chandrashekar(up::AbstractArray, un::AbstractArray, param::parameters) # copied from (Chan 2018)
    if param.dim == 1
        f1 = zeros(size(up))

        @views velp = up[:,2:param.dim+1] ./ up[:,1]
        @views veln = un[:,2:param.dim+1] ./ un[:,1]
        velavg = 0.5 .* (velp .+ veln)

        @views betap = 0.5 .* up[:,1] ./ Euler_pressure(up, param)
        @views betan = 0.5 .* un[:,1] ./ Euler_pressure(un, param)

        @views f1[:,1] .= logmean.(up[:,1], un[:,1]) .* velavg[:,1]
        @views f1[:,2] .= 0.5 .* (up[:,1] .+ un[:,1]) ./ (betap .+ betan) .+ velavg .* f1[:,1]
        @views f1[:,3] .= f1[:,1] .* (0.5 / (param.gamma-1) ./ logmean.(betan,betap) .- 0.25 .* (velp.^2 .+ veln.^2)) .+ velavg .* f1[:,2]

        return [f1]
    else
        error("Invalid number of dimensions.")
    end
end
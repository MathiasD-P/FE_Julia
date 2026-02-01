#####################################################################
# Everything concerning the problem's physics
#####################################################################

# THE SOLUTION IS ALWAYS A MATRIX OF SIZE (Nstate, NDOF)

#####################################################################
# Physical Fluxes
#####################################################################

function compute_physflux!(f, u, param::parameters)
    if param.pdetype == "LinAdv"
        @. f[1] = param.a * u
        return f

    elseif param.pdetype == "Burgers"
        @. f[1] = 0.5 * u^2
        return f

    elseif param.pdetype == "EulerPerfGas"
        for index in axes(u,1)
            @views localflux = Euler_physflux(u[index,:], param)
            for dir in 1:param.dim
                f[dir][index,:] .= localflux[dir]
            end
        end

        return f
    end
end

#####################################################################
# Numerical Fluxes
#####################################################################

function compute_numflux(f, un, up, nphys, param::parameters)
    if param.pdetype == "LinAdv"
        if param.numfluxtype == "central"
            @. f[1] = 0.5 * param.a * (up + un)
            return f
        elseif param.numfluxtype == "upwind"
            f[1] .= up
            if param.a * nphys[1] >= 0
                index = nphys .== nphys[1]
            else
                index = nphys .!= nphys[1]
            end
            f[1][index] = un[index]
            @. f[1] = param.a * f[1]
            return f
        else
            print("Undefined numerical flux!")
        end

    elseif param.pdetype == "Burgers"
        if param.numfluxtype == "central"
            @. f[1] = 0.25 * (up^2 + un^2)
            return f
        elseif param.numfluxtype == "LF"
            @. f[1] = 0.25 * ((un^2 + up^2) - max(abs(up), abs(un)) * (up - un) * dg.nphys)
            return f
        elseif param.numfluxtype == "EC_split"
            @. f[1] = (1/6) * (un^2 + up * un + up^2)
            return f
        else
            print("Undefined numerical flux!")
        end

    elseif param.pdetype == "EulerPerfGas"        
        for index in axes(up,1)
            if param.numfluxtype == "central"
                @views localflux = @. 0.5*(Euler_physflux(un[index,:], param) + Euler_physflux(up[index,:], param))
            elseif param.numfluxtype == "EC_Chandrashekar"
                @views localflux = Euler_numflux_Chandrashekar(up[index,:], un[index,:], param)
            else
                print("Undefined numerical flux!")
            end

            for dir in 1:param.dim
                f[dir][index,:] .= localflux[dir]
            end
        end

        return f
    end
end

#####################################################################
# Two-point Fluxes
#####################################################################

function two_pt_flux(u, uf, param::parameters)
    M, N = size(u)
    Npts = M + size(uf, 1)
    nphys = nothing # We don't need the physical normal here, but this is a required argument for compute_numflux

    if param.dim == 1
        F1 = zeros((Npts,Npts,N))
        for j in 1:Npts
            for i in 1:j
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

                f = compute_numflux(up,un, nphys, param)
                F1[i,j,:] .= f[1][:]
                F1[j,i,:] .= f[1][:]
            end
        end

        return [F1]

    elseif param.dim == 2
        F1 = zeros((Npts,Npts,N))
        F2 = zeros((Npts,Npts,N))
        for j in 1:Npts
            for i in 1:j
                if j > M
                    up = @view uf[j - M,:]
                else
                    up = @view u[j,:]
                end
                if i > M
                    un = @view uf[i - M,:]
                else
                    un = @view u[i,:]
                end

                f = compute_numflux(up,un, nphys, param)
                F1[i,j,:] .= f[1][:]
                F1[j,i,:] .= f[1][:]
                F2[i,j,:] .= f[2][:]
                F2[j,i,:] .= f[2][:]
            end
        end

        return [F1, F2]
    end
end

#####################################################################
# Entropy mappings
#####################################################################

function compute_evar(u, param::parameters)
    if param.pdetype == "Burgers"
        return u
    elseif param.pdetype == "EulerPerfGas"
        v = zeros(size(u))
        for index in axes(u,1)
            @views v[index,:] .= Euler_evar(u[index,:], param)
        end

        return v
    end
end

function compute_cvar(v, param::parameters)
    if param.pdetype == "Burgers"
        return v
    elseif param.pdetype == "EulerPerfGas"
        u = zeros(size(u))
        for index in axes(u,1)
            @views u[index,:] .= Euler_evar(u[index,:], param)
        end

        return u
    end
end

#####################################################################
# Compute Entropy
#####################################################################

function compute_local_entropy(u, param::parameters)
    if param.pdetype == "LinAdv"
        return 0.5 .* u.^2
    elseif param.pdetype == "Burgers"
        return 0.5 .* u.^2
    elseif param.pdetype == "EulerPerfGas"
        s = zeros((size(u,1),1))
        for index in axes(u,1)
            @views s[index] .= -u[index,1] .* Euler_cvar_entropy(u[index,:], param)
        end

        return s
    end
end

#####################################################################
# Euler Helper functions
# These apply to one SINGLE NODE
#####################################################################

# (\rho * e)(u)
function Euler_cvar_intenergy(u, param::parameters)
    if param.dim == 1
        return u[end] - 0.5 * u[2]^2 / u[1]
    elseif param.dim == 2
        return u[end] - 0.5 * (u[2]^2 + u[3]^2) / u[1]
    end
end

# (\rho * e)(v)
function Euler_evar_intenergy(v, param::parameters)
    return ((param.gamma-1) / (-v[end])^param.gamma)^(1/(param.gamma-1)) * exp(-Euler_evar_entropy(v,param) / (param.gamma-1))
end

# (s / cv)(u)
function Euler_cvar_entropy(u, param::parameters)
    return log(Euler_pressure(u, param) / u[1]^param.gamma)
end

# (s / cv)(v)
function Euler_evar_entropy(v, param::parameters)
    if param.dim == 1
        return param.gamma - v[1] + 0.5 * v[2]^2 / v[end]
    elseif param.dim == 2
        return param.gamma - v[1] + 0.5 * (v[2]^2 + v[3]^2) ./ v[end]
    end
end

# p(u)
function Euler_pressure(u, param::parameters)
    return (param.gamma - 1) * Euler_cvar_intenergy(u, param)
end

# v(u)
function Euler_evar(u, param::parameters)
    rhoe = Euler_cvar_intenergy(u,param)
    s = Euler_cvar_entropy(u, param)

    v = zeros(size(u))
    v[1] = (-s + param.gamma + 1) - u[end] / rhoe
    v[2:1+param.dim] .= u[2:1+param.dim] ./ rhoe
    v[end] = -u[1] / rhoe

    return v
end

# u(v)
function Euler_cvar(v, param::parameters)
    rhoe = Euler_evar_intenergy(v, param)

    u = zeros(size(v))
    u[1] = -rhoe * v[end]
    u[2:1+param.dim] .= v[2:1+param.dim] .* rhoe
    u[end] = rhoe * (1 - 0.5 * sum(v[:,2:1+param.dim].^2) /  v[end])

    return u
end

# f_dir
function Euler_physflux(u, param::parameters)
    p = Euler_pressure(u, param)

    f1 = zeros(size(u))
    f1[1] = u[2]
    @views f1[2:param.dim+1] .=  u[2:param.dim+1] .* u[2] ./  u[1]
    f1[2] = f1[2] + p
    f1[end] = (p +  u[end]) * u[2] / u[1]

    if param.dim == 1
        return [f1]
    elseif param.dim == 2
        f2 = zeros(size(u))
        f2[1] = u[3]
        f2[2] = f1[3]
        f2[3] = u[3]^2 /  u[1] + p
        f2[end] = (p +  u[end]) * (u[3]/ u[1])

        return [f1, f2]
    end

end

# Numerical fluxes for Euler

function logmean(up, un)
    epsilon = 1e-2 # tolerance for logmean computation

    if abs(up - un) < epsilon # Roe simplification for up -> un
        r = up/un
        a = ((r - 1) / (r + 1))^2
        return 0.5 * (up + un) / (1 + a/3 + a^2/5 + a^3/7)
    else
        return (up - un) / (log(up) - log(un))
    end
end

# ONLY FOR 1D RIGHT NOW
function Euler_numflux_Chandrashekar(up, un, param::parameters) # copied from (Chan 2018)
    if param.dim == 1
        f1 = zeros(size(up))

        velp = up[2] / up[1]
        veln = un[2] / un[1]
        velavg = 0.5 * (velp + veln)

        betap = 0.5 * up[1] / Euler_pressure(up, param)
        betan = 0.5 * un[1] / Euler_pressure(un, param)

        f1[1] = logmean(up[1], un[1]) * velavg[1]
        f1[2] = 0.5 * (up[1] + un[1]) / (betap + betan) + velavg * f1[1]
        f1[3] = f1[1] * (0.5 / (param.gamma-1) / logmean(betan,betap) - 0.25 * (velp^2 + veln^2)) + velavg[1] * f1[2]

        return [f1]
    else
        error("Invalid number of dimensions.")
    end
end
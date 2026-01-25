#####################################################################
# Everything concerning the problem's physics
#####################################################################

using LinearAlgebra

# THE SOLUTION IS ALWAYS A MATRIX OF SIZE (Nstate, NDOF)

#####################################################################
# Physical Fluxes
#####################################################################

function compute_physflux(u, param::parameters)
    if param.pdetype == "LinAdv"
        return [param.a .* u]

    elseif param.pdetype == "Burgers"
        return [0.5 .* u.^2]

    elseif param.pdetype == "EulerPerfGas"
        return Euler_physflux(u, param)
    end
end

#####################################################################
# Numerical Fluxes
#####################################################################

function compute_numflux(un, up, nphys, param::parameters)
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
            print("Undefined numerical flux!")
        end

    elseif param.pdetype == "Burgers"
        if param.numfluxtype == "central"
            return [0.25 .* (up.^2 .+ un.^2)]
        elseif param.numfluxtype == "LF"
            return [0.25 .* ((un.^2 .+ up.^2) .- max.(abs.(up), abs.(un)) .* (up .- un) .* dg.nphys)]
        elseif param.numfluxtype == "EC_split"
            return [(1/6) .* (un.^2 .+ up .* un .+ up.^2)]
        else
            print("Undefined numerical flux!")
        end

    elseif param.pdetype == "EulerPerfGas"
        if param.numfluxtype == "central"
            return 0.5 .* (Euler_physflux(up, param) .+ Euler_physflux(un, param))
        elseif param.numfluxtype == "EC_Chandrashekar"
            return Euler_numflux_Chandrashekar(up, un, param)
        else
            print("Undefined numerical flux!")
        end
    end
end

#####################################################################
# Two-point Fluxes
#####################################################################

function two_pt_flux(up, un, param::parameters)
    M,N = size(up)

    if param.dim == 1
        F1 = zeros((N,N,M))
        for j in axes(up, 1)
            for i in 1:j
                upv = @view up[i,:]
                unv = @view un[j,:]
                f = numflux(upv,unv,param)
                F1[i,j,:] .= f[1]
                F1[i,j,:] .= f[1]
            end
        end

        return [F1]

    elseif param.dim == 2
        F1 = zeros((N,N,M))
        F2 = zeros((N,N,M))
        for j in axes(up, 1)
            for i in 1:j
                upv = @view up[i,:]
                unv = @view un[j,:]
                f = numflux(upv,unv,param)
                F1[i,j,:] .= f[1]
                F1[j,i,:] .= f[1]
                F2[i,j,:] .= f[2]
                F2[j,i,:] .= f[2]
            end
        end

        return [F1, F2]
    end
end

#####################################################################
# Entropy mappings
#####################################################################

function eval_evar(u, param::parameters)
    if param.pdetype == "Burgers"
        return u
    elseif param.pdetype == "EulerPerfGas"
        rhoe = cvar_intenergy(u,param)
        s = cvar_entropy(u, param)

        v = zeros(size(u))
        v[1,:] .= (-s .+ param.gamma .+ 1) .- u[end,:] ./ rhoe
        v[2:1+param.dim,:] .= u[2:1+param.dim,:] ./ rhoe
        v[end,:] .= -u[1,:] ./ rhoe

        return v
    end
end

function eval_cvar(v, param::parameters)
    if param.pdetype == "Burgers"
        return v
    elseif param.pdetype == "EulerPerfGas"
        rhoe = evar_intenergy(v, param)

        u = zeros(size(u))
        u[:,1] .= -rhoe .* v[:,end]
        u[:,2:1+param.dim] .= v[:,2:1+param.dim] .* rhoe
        u[:,end] .= -v[:,1] .+ param.gamma .+ 0.5 .* sum(v[:,2:1+param.dim].^2, dims=2)' ./  v[:,end]

        return u
    end
end

#####################################################################
# Euler helper functions
#####################################################################

# (\rho * e)(u)
function cvar_intenergy(u, param::parameters)
    if param.dim == 1
        return u[:,end] .- 0.5 .* u[:,2].^2 ./ u[:,1]
    elseif param.dim == 2
        return u[:,end] .- 0.5 .* (u[:,2].^2 .+ u[3,:].^2) ./ u[1,:]
    end
end

# (\rho * e)(v)
function evar_intenergy(v, param::parameters)
    return ((param.gamma-1) ./ (-v[:,end]).^param.gamma).^(1/(param.gamma-1)) .* exp.(-evar_entropy(v,param) ./ (param.gamma-1))
end

# (s / cv)(u)
function cvar_entropy(u, param::parameters)
    return log.(pressure(u, param) ./ u[:,1].^param.gamma)
end

# (s / cv)(v)
function evar_entropy(v, param::parameters)
    if param.dim == 1
        return param.gamma .- v[:,1] .+ 0.5 .* v[:,2].^2 ./ v[:,end]
    elseif param.dim == 2
        return param.gamma .- v[:,1] .+ 0.5 .* (v[:,2].^2 .+ v[:,3].^2) ./ v[:,end]
    end
end

# p(u)
function pressure(u, param::parameters)
    return (param.gamma - 1) .* cvar_intenergy(u, param)
end

# logmean
function logmean(up, un)
    epsilon = 1e-2

    if abs(up - un) < epsilon # Roe simplification for up -> un
        a = ((up/un - 1) / (up/un + 1))^2
        return 0.5 * (up + un) / (1 + a/3 + a^2/5 + u^3/7)
    else
        return (up - un) / (log(up) - log(un))
    end
end

# f_dir
function Euler_physflux(u, param::parameters)
    p = pressure(u, param)

    f1 = zeros(size(u))
    f1[:,1] .= u[:,2]
    f1[:2:param.dim+1] .=  u[:,2:param.dim+1] .* u[:,2]' ./  u[:,1]
    f1[:,2] .= f[:,2] .+ p
    f1[:,end] .= (p .+  u[:,end]) .* (u[:,2]./ u[:,1])

    if param.dim == 1
        return [f1]
    elseif param.dim == 2
        f2 = zeros(size(u))
        f2[:,1] .= u[3,:]
        f2[:,2] .= f1[:,3]
        f2[:,3] .= u[:,3].^2 ./  u[:,1] .+ p
        f2[:,end] .= (p .+  u[:,end]) .* (u[:,3]./ u[:,1])

        return [f1, f2]
    end

end

function Euler_numflux_Chandrashekar(up, un, param::parameters)
    uavg2 = 0.5 .* (sum((up[:,2:param.dim+1] .+ un[:,2:param.dim+1]).^2, dims=2) .- sum(up[:,2:param.dim+1].^2 .+ un[:,2:param.dim+1].^2, dims=2))
    uavg = 0.5 .* (up[:,2:param.dim+1] ./ up[:,1] .+ un[:,2:param.dim+1] ./ un[:,1])
    pl = pressure(up, param)
    pr = pressure(un, param)
    rholog = logmean.(up[1,:], un[1,:])
    pavg = (up[:,1].+un[:,1]) ./ (up[:,1] ./ pl .+ un[:,1] ./ pr)
    plog = 0.5 .* rholog ./ logmean.(0.5 .* up[:,1] ./ pl, 0.5 .* un[:,1] ./ pr)

    f1 = zeros(size(up))
    f1[:,1] .= rholog .* uavg[:,1]
    f1[:,2:param.dim+1] .= uavg .* f1[:,1]
    f1[:,2] .= f1[:,2] .+ pavg
    f1[:,end] .= (plog ./ (param.gamma-1) .+ pavg .+ 0.5 .* uavg2) .* uavg[:,1]

    if param.dim == 1
        return [f1]
    elseif param.dim == 2
        f2 = zeros(size(up))
        f2[:,1] .= rholog .* uavg[:,2]
        f2[:,2] .= f1[:,3]
        f2[:,3] .= f2[:,1] .* uavg[:,2] .+ pavg
        f1[:,end] .= (plog ./ (param.gamma-1) .+ pavg .+ 0.5 .* uavg2) .* uavg[:,2]

        return [f1, f2]
    end
end
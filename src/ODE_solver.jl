#####################################################################
# Performs explicit time-stepping
#####################################################################

function ODE_solver(u0::Matrix{Float64}, BChandler::Dict, dg::DG, param::parameters)

    if param.calc_entropy
        S = zeros((param.Nsteps,)) # pre-allocate memory for entropy calculation if required.
    end

    if param.ODE_solver == "LSERK45"
        RKa = ( 0.0,
                -567301805773.0/1357537059087.0,
                -2404267990393.0/2016746695238.0,
                -3550918686646.0/2091501179385.0,
                -1275806237668.0/842570457699.0);
        RKb = ( 1432997174477.0/9575080441755.0,
                5161836677717.0/13612068292357.0,
                1720146321549.0/2090206949498.0,
                3134564353537.0/4481467310338.0,
                2277821191437.0/14882151754819.0);
        RKc = ( 0.0,
                1432997174477.0/9575080441755.0,
                2526269341429.0/6820363962896.0,
                2006345519317.0/3224310063776.0,
                2802321613138.0/2924317926251.0);
        RKstages = 5

        u = u0
        current_time = 0
        residual = zeros(size(u))
        rhs = zeros(size(u))

        for istep in 1:param.Nsteps
            for stage in 1:RKstages

                rktime = current_time + RKc[stage] * param.dt

                rhs = build_residual!(rhs, u, rktime, BChandler, dg, param)

                residual .= RKa[stage] .* residual .+ param.dt .* rhs
                u .= u .+ RKb[stage] .* residual
            end
            current_time += param.dt

            if param.calc_entropy
                S[istep] = compute_total_entropy(u, dg, param)
            end
        end

        if param.calc_entropy
            return (u, current_time, S)
        else
            return (u, current_time)
        end

    end
end

function compute_total_entropy(u, dg::DG, param::parameters)
    if dg.mesh isa LMesh
        s = compute_local_entropy(block_matmul(dg.refelem.chiq, u, dg.mesh.Nel), param)
        s = block_matmul(Diagonal(dg.refelem.wq), s, dg.mesh.Nel)

        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nqnodes*(ielem-1):dg.refelem.Nqnodes*ielem
            @views s[index,:] .= s[index,:] .* dg.mesh.J[ielem]
        end

        return sum(s)
    end
end
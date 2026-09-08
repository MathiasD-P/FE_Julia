#####################################################################
# Performs explicit time-stepping
#####################################################################

abstract type TimeIntegrator end

struct AuxiliaryCalc
    calc_entropy::Bool
end

struct LSERK45 <: TimeIntegrator
    # Time Integration settings
    Nsteps::Int
    dt::Float64

    # Auxiliary calculations
    auxcalc::AuxiliaryCalc
end

struct SSPRK33 <: TimeIntegrator
    # Time Integration settings
    Nsteps::Int
    dt::Float64

    # Auxiliary calculations
    auxcalc::AuxiliaryCalc
end

struct SSPRRK33 <: TimeIntegrator
    # Time Integration settings
    Nsteps::Int
    dt::Float64

    # Auxiliary calculations
    auxcalc::AuxiliaryCalc
end

function ODE_solver(u0::Matrix{Float64}, dg::DG, physics::PhysProp, method::LSERK45)
    if method.auxcalc.calc_entropy
        S = zeros((method.Nsteps,)) # pre-allocate memory for entropy calculation if required.
    end

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

    # Required initializations for LSERK45
    u = u0
    current_time = 0.0
    residual = zeros(size(u))
    rhs = zeros(size(u))

    for istep in 1:method.Nsteps
        for stage in 1:RKstages

            rktime = current_time + RKc[stage] * method.dt

            rhs = build_residual!(rhs, u, rktime, dg, physics)

            residual .= RKa[stage] .* residual .+ method.dt .* rhs
            u .= u .+ RKb[stage] .* residual
        end
        current_time += method.dt

        if method.auxcalc.calc_entropy
            S[istep] = compute_total_entropy(u, dg, physics.PDE)
        end
    end

    if method.auxcalc.calc_entropy
        return Dict("solution" => u, "time" => current_time, "entropy" => S)
    else
        return Dict("solution" => u, "time" => current_time)
    end
end

function ODE_solver(u0::Matrix{Float64}, dg::DG, physics::PhysProp, method::SSPRK33)
    if method.auxcalc.calc_entropy
        S = zeros((method.Nsteps,)) # pre-allocate memory for entropy calculation if required.
    end

    # Required initializations for SSPRK33
    ubuff = copy(u0)
    u = u0
    current_time = 0.0
    rhs = Matrix{Float64}(undef, size(u0)...)

    for istep in 1:method.Nsteps
        rhs = build_residual!(rhs, u, current_time, dg, physics)
        @. ubuff = u + method.dt * rhs

        rhs = build_residual!(rhs, ubuff, current_time + method.dt, dg, physics)
        @. ubuff = 0.75 * u + 0.25 * ubuff + 0.25 * method.dt * rhs

        rhs = build_residual!(rhs, ubuff, current_time + 0.5 * method.dt, dg, physics)
        @. u = (1/3) * u + (2/3) * ubuff + (2/3) * method.dt * rhs

        if method.auxcalc.calc_entropy
            S[istep] = compute_total_entropy(u, dg, physics.PDE)
        end

        current_time += method.dt
    end

    if method.auxcalc.calc_entropy
        return Dict("solution" => u, "time" => current_time, "entropy" => S)
    else
        return Dict("solution" => u, "time" => current_time)
    end
end

function ODE_solver(u0::Matrix{Float64}, dg::DG, physics::PhysProp, method::SSPRRK33)
    if method.auxcalc.calc_entropy
        S = zeros((method.Nsteps,)) # pre-allocate memory for entropy calculation if required.
    end

    if !(physics.PDE isa Burgers)
        error("Relaxation RK is only implemented for Burgers at the moment!")
    end

    RKa = [0.0  0.0;
            1.0  0.0;
            0.25 0.25];     
    RKb = (1/6,
            1/6,
            2/3);
    RKc = (0.0,
            1.0,
            0.5);
    RKstages = 3

    # Required initializations for SSPRK33
    rhs = Tuple(Matrix{Float64}(undef, size(u0)...) for s in 1:3)
    u = u0
    ubuff = copy(u0)
    current_time = 0.0

    for istep in 1:method.Nsteps
        build_residual!(rhs[1], u, current_time, dg, physics)
        gammanum = 0.0
        gammaden = 0.0

        for stage in 2:RKstages
            ubuff .= u
            for istage in 1:stage-1
                @. ubuff += method.dt * RKa[stage, istage] * rhs[istage]
            end

            build_residual!(rhs[stage], ubuff, current_time + method.dt * RKc[stage], dg, physics)
        end

        for istage in 1:RKstages
            for jstage in 1:istage
                if istage == jstage
                    gammaden += RKb[istage] * RKb[jstage] * dot(rhs[istage], block_matmul(dg.refelem.M, rhs[jstage], dg.mesh.Nel))
                else
                    delta = dot(rhs[istage], block_matmul(dg.refelem.M, rhs[jstage], dg.mesh.Nel))
                    gammaden += 2 * RKb[istage] * RKb[jstage] * delta
                    gammanum += RKb[istage] * RKa[istage, jstage] * delta
                end
            end
        end

        if abs(gammaden) < 1e-15
            gamma = 1.0
        else
            gamma = 2.0 * gammanum / gammaden
        end

        for stage in 1:RKstages
            @. u += gamma * method.dt * RKb[stage] * rhs[stage]
        end

        current_time += method.dt

        if method.auxcalc.calc_entropy
            S[istep] = compute_total_entropy(u, dg, physics.PDE)
        end
    end

    if method.auxcalc.calc_entropy
        return Dict("solution" => u, "time" => current_time, "entropy" => S)
    else
        return Dict("solution" => u, "time" => current_time)
    end
end

function compute_total_entropy(u, dg::DG, PDE::GoverningPDE)
    if dg.mesh isa LMesh
        s = compute_local_entropy(block_matmul(dg.refelem.chiq, u, dg.mesh.Nel), PDE)
        s = block_matmul(Diagonal(dg.refelem.wq), s, dg.mesh.Nel)

        for ielem = 1:dg.mesh.Nel
            index = 1+dg.refelem.Nqnodes*(ielem-1):dg.refelem.Nqnodes*ielem
            @views s[index,:] .= s[index,:] .* dg.mesh.J[ielem]
        end

        return sum(s)
    end
end

# function ODE_solver(u0::Matrix{Float64}, BChandler::Dict, dg::DG, param::parameters)

#     if param.calc_entropy
#         S = zeros((param.Nsteps,)) # pre-allocate memory for entropy calculation if required.
#     end

#     if param.ODE_solver == "LSERK45"
#         RKa = ( 0.0,
#                 -567301805773.0/1357537059087.0,
#                 -2404267990393.0/2016746695238.0,
#                 -3550918686646.0/2091501179385.0,
#                 -1275806237668.0/842570457699.0);
#         RKb = ( 1432997174477.0/9575080441755.0,
#                 5161836677717.0/13612068292357.0,
#                 1720146321549.0/2090206949498.0,
#                 3134564353537.0/4481467310338.0,
#                 2277821191437.0/14882151754819.0);
#         RKc = ( 0.0,
#                 1432997174477.0/9575080441755.0,
#                 2526269341429.0/6820363962896.0,
#                 2006345519317.0/3224310063776.0,
#                 2802321613138.0/2924317926251.0);
#         RKstages = 5

#         # Required initializations for LSERK45
#         u = u0
#         current_time = 0.0
#         residual = zeros(size(u))
#         rhs = zeros(size(u))

#         for istep in 1:param.Nsteps
#             for stage in 1:RKstages

#                 rktime = current_time + RKc[stage] * param.dt

#                 rhs = build_residual!(rhs, u, rktime, BChandler, dg, param)

#                 residual .= RKa[stage] .* residual .+ param.dt .* rhs
#                 u .= u .+ RKb[stage] .* residual
#             end
#             current_time += param.dt

#             if param.calc_entropy
#                 S[istep] = compute_total_entropy(u, dg, param)
#             end
#         end

#     elseif param.ODE_solver == "SSPRK33"
#         # Required initializations for SSPRK33
#         ubuff = copy(u0)
#         u = u0
#         current_time = 0.0
#         rhs = Matrix{Float64}(undef, size(u0)...)

#         for istep in 1:param.Nsteps
#             rhs = build_residual!(rhs, u, current_time, BChandler, dg, param)
#             @. ubuff = u + param.dt * rhs

#             rhs = build_residual!(rhs, ubuff, current_time + param.dt, BChandler, dg, param)
#             @. ubuff = 0.75 * u + 0.25 * ubuff + 0.25 * param.dt * rhs

#             rhs = build_residual!(rhs, ubuff, current_time + 0.5 * param.dt, BChandler, dg, param)
#             @. u = (1/3) * u + (2/3) * ubuff + (2/3) * param.dt * rhs

#             if param.calc_entropy
#                 S[istep] = compute_total_entropy(u, dg, param)
#             end

#             current_time += param.dt
#         end

#     elseif param.ODE_solver == "SSPRRK33" #NOT EFFICIENT JUST FOR TESTING, TEMPORAL DICRETIZATION IS NOT THE FOCUS OF OUR PROJECT
#         if param.pdetype != "Burgers"
#             error("Relaxation RK is only implemented for Burgers at the moment!")
#         end

#         RKa = [0.0  0.0;
#                1.0  0.0;
#                0.25 0.25];     
#         RKb = (1/6,
#                1/6,
#                2/3);
#         RKc = (0.0,
#                1.0,
#                0.5);
#         RKstages = 3

#         # Required initializations for SSPRK33
#         rhs = Tuple(Matrix{Float64}(undef, size(u0)...) for s in 1:3)
#         u = u0
#         ubuff = copy(u0)
#         current_time = 0.0

#         for istep in 1:param.Nsteps
#             build_residual!(rhs[1], u, current_time, BChandler, dg, param)
#             gammanum = 0.0
#             gammaden = 0.0

#             for stage in 2:RKstages
#                 ubuff .= u
#                 for istage in 1:stage-1
#                     @. ubuff += param.dt * RKa[stage, istage] * rhs[istage]
#                 end

#                 build_residual!(rhs[stage], ubuff, current_time + param.dt * RKc[stage], BChandler, dg, param)
#             end

#             for istage in 1:RKstages
#                 for jstage in 1:istage
#                     if istage == jstage
#                         gammaden += RKb[istage] * RKb[jstage] * dot(rhs[istage], block_matmul(dg.refelem.M, rhs[jstage], dg.mesh.Nel))
#                     else
#                         delta = dot(rhs[istage], block_matmul(dg.refelem.M, rhs[jstage], dg.mesh.Nel))
#                         gammaden += 2 * RKb[istage] * RKb[jstage] * delta
#                         gammanum += RKb[istage] * RKa[istage, jstage] * delta
#                     end
#                 end
#             end

#             if abs(gammaden) < 1e-15
#                 gamma = 1.0
#             else
#                 gamma = 2.0 * gammanum / gammaden
#             end

#             for stage in 1:RKstages
#                 @. u += gamma * param.dt * RKb[stage] * rhs[stage]
#             end

#             current_time += param.dt

#             if param.calc_entropy
#                 S[istep] = compute_total_entropy(u, dg, param)
#             end
#         end
#     end

#     if param.calc_entropy
#         return Dict("solution" => u, "time" => current_time, "entropy" => S)
#     else
#         return Dict("solution" => u, "time" => current_time)
#     end
# end
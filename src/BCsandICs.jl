function initialize_states(dg::DG, param::parameters, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    u = zeros((size(pts,1), dg.Nstates))

    if param.ICname == "sin_1state"
        u .= sin.(2*pi .* pts)

        return u
    
    elseif param.ICname == "exp_1state"
        u .= exp.(-80 .* pts.^2)

        return u
    
    elseif param.ICname == "IsentropicDensityWave" # Classical test case used by (Chan 2025)
        A = 0.5 # tunable wave amplitude
        u[:,1] .= 1 .+ A .* sin.(2*pi.*sum(pts, dims=2))
        @views u[:,2] .= 0.1 .* u[:,1]

        if param.dim == 1
            @views u[:,end] = 10 / (param.gamma - 1) .+ (0.5 * 0.1^2) ./ u[:,1]
        elseif param.dim == 2
            @views u[:,] .= 0.2 .* u[:,1]
            @views u[:,end] = 10 / (param.gamma - 1) .+ (0.5 * (0.1^2 +0.2^2)) ./ u[:,1]
        end

        return u


    elseif param.ICname == "SodShockTube"
        u[pts .< 0,1] .= 1
        u[pts .>= 0,1] .= 0.125
        u[pts .< 0,3] .= 1.0 / (param.gamma-1)
        u[pts .< 0,3] .= 0.1 / (param.gamma-1)

        return u
    end
end

# consistent mesh and BCHandler initializations (must absolutely be updated together)

function initialize_mesh(param::parameters)
    if param.domain == "unit_square_linear_quad"
        if param.BCname == "periodic"
            return make_rectangle_quad(param.Neldim, param.Neldim)
        elseif param.BCname == "homogeneous_Dirichlet"
            return make_rectangle_quad(Neldim, Neldim, [-1,-1,-1,-1])
        end

    elseif param.domain == "unit_interval_linear"
        if param.BCname == "periodic"
            return make_interval(collect(range(-0.5, 0.5, param.Neldim+1)), [0, 0])
        elseif param.BCname == "homogeneous_Dirichlet"
            return make_interval(collect(range(-0.5, 0.5, param.Neldim+1)), [-1, -1])
        elseif param.BCname == "SodShockTube"
            return make_interval(collect(range(-0.5, 0.5, param.Neldim+1)), [-1, -2])
        end
    end
end

function initialize_BCHandler(dg::DG, param::parameters)
    if param.BCname == "periodic"
        return Dict()

    elseif param.BCname == "homogeneous_Dirichlet"
        return Dict(-1 => zeros((dg.Nstates,)))

    elseif param.BCname == "SodShockTube"
        return Dict(-1 => [1.0, 0.0, 1.0 / (param.gamma-1)], -2 => [0.125, 0.0, 0.1 / (param.gamma-1)])
    end
end

function evaluate_BC(BCHandler::Dict, dg::DG, t)
    if isempty(BCHandler)
        return spzeros(dg.NFval, dg.Nstates)

    elseif BCHandler isa Dict{<:Integer, <:Vector{<:Real}}
        return sum([dg.BFtoF[itag] * BCHandler[dg.mesh.BCtags[itag]]' for itag in 1:dg.mesh.Ntags])
    end
end
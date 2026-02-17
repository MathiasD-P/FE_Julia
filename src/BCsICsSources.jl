function initialize_states(dg::DG, param::parameters, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    u = zeros((size(pts,1), dg.Nstates))

    if param.ICname == "sin_1state"
        if param.dim != 1
            error("IC only works in 1D!")
        end

        u .= sin.(2*pi .* pts)

        return u
    
    elseif param.ICname == "exp_1state"
        if param.dim != 1
            error("IC only works in 1D!")
        end

        u .= exp.(-80 .* pts.^2)

        return u
    
    elseif param.ICname == "IsentropicDensityWave" # Classical test case used by (Chan 2025)
        A = 0.5 # tunable wave amplitude
        u[:,1] .= 1 .+ A .* sin.(2*pi.*sum(pts, dims=2))
        @views u[:,2] .= 0.1 .* u[:,1]

        if param.dim == 1
            @views u[:,end] = 10 / (param.gamma - 1) .+ (0.5 * 0.1^2) .* u[:,1]
        elseif param.dim == 2
            @views u[:,3] .= 0.2 .* u[:,1]
            @views u[:,end] = 10 / (param.gamma - 1) .+ (0.5 * (0.1^2 +0.2^2)) .* u[:,1]
        end

        return u
    
    elseif param.ICname == "GaussianVelocity"
        u[:,1] .= 1.0
        u[:,2:end-1] .= exp.(-80 .* sum(pts.^2, dims=2)) .+ 1.0
        u[:,3] .= 5.0

        return u

    elseif param.ICname == "SodShockTube"
        if param.dim != 1
            error("Sod shock only works in 1D!")
        end

        u[pts .< 0,1] .= 1
        u[pts .>= 0,1] .= 0.125
        u[pts .< 0,3] .= 1.0 / (param.gamma-1)
        u[pts .< 0,3] .= 0.1 / (param.gamma-1)

        return u
    
    # MANUFACTURED SOLUTIONS
    elseif param.ICname != param.sourcename
        error("Source and manufactured initial conditions must match!")
    
    elseif param.ICname == "GassnerEuler"
        u[:,1] .= 2 .+ 0.1 .* sin.(2*pi .* sum(pts, dims=2))
        @views u[:,2] .= u[:,1]
        @views u[:,3] .= u[:,1].^2

        return u
    
    else
        error("Unknown IC name!")
    end
end

function compute_source(dg::DG, param::parameters, pts::Matrix{Float64}, time::Float64)
    Q = Matrix{Float64}(undef, (size(pts,1), dg.Nstates))

    if param.dim == 2
        println("Verifiy the accuracy of source and manufactured solution in 2D. This has not been validated...")
    end

    if param.sourcename == "GassnerEuler" # not exactly the same as in Gassner, corrected source and solution from CPerthick (should work for 1 and 2D)
        Q[:,1] .= - 2 * pi * 0.1 * cos.(2*pi * (sum(pts, dims=2) .- 2*time))
        Q[:,2:end-1] .= 1/100 * 2 * pi .* cos.(2*pi .* (sum(pts, dims=2) .- 2*time)) .* (35 * param.gamma .- 45 .+ 2 .* (param.gamma - 1) .* sin.(2*pi .* (sum(pts, dims=2) .- 2*time)))
        Q[:,end] .= 1/100 * 2 * pi * cos.(2*pi * (sum(pts, dims=2) .- 2*time)) .* (-75 .+ 35 * param.gamma .+ 2 * (param.gamma - 2) .* sin.(2*pi .* (sum(pts, dims=2) .- 2*time)))
    end

    return Q
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
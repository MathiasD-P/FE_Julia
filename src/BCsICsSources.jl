function initialize_states(dg::DG, param::parameters, pts = nothing)
    flag_bpts = false

    if isnothing(pts)
        pts = dg.bpts
        flag_bpts = true
    end

    u = zeros((size(pts,1), dg.Nstates))

    if param.ICname == "sin_1state"
        if param.dim != 1
            error("IC only works in 1D!")
        end

        u .= sin.(2*pi .* param.k[1] .* (pts .- param.phi)) .+ param.av

        return u
    
    elseif param.ICname == "exp_1state"
        if param.dim != 1
            error("IC only works in 1D!")
        end

        u .= exp.(-80 .* pts.^2)

        return u
    
    elseif param.ICname == "rarefaction_1state"
        if param.dim != 1
            error("IC only works in 1D!")
        end

        u = ones(size(pts))
        u[pts .<= 0.0] .= -1.0

        return u
    
    elseif param.ICname == "zeros_and_ones"
        if param.dim != 1
            error("IC only works in 1D!")
        end

        u = [Float64(isodd(i)) for i in 1:length(pts)] 
        u = reshape(u, (length(u),1))

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
    
    elseif param.ICname == "ChanWave" # Scaled test case used by (Chan 2025) to investigate convergence of viscosity and entropy deficit.
        if param.dim != 1
            error("Only works for 1D right now!")
        end

        A = 0.5 # tunable wave amplitude
        u[:,1] .= 1 .+ A .* sin.(2*pi.*pts .+ 0.1/2)
        u[:,2] .= A .* sin.(2*pi.*pts .+ 0.2/2) .* u[:,1]
        u[:,3] .= u[:,1].^param.gamma ./ (param.gamma-1) .+ 0.5 .* u[:,2].^2 ./ u[:,1]

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
    
    elseif param.ICname == "Burgulence"
        if !flag_bpts
            error("Can only initialize Burgulence at basis nodes!")
        end

        if (param.domain != "unit_interval_linear") || (param.BCname != "periodic")
            error("Burgulence can only be initialized on the unit circle!")
        end

        bnodes = make_nodes(param.bnodes)
        Burg = Burgulence(bnodes, param.Neldim, param.kmax)

        return Burg.IC

    
    # MANUFACTURED SOLUTIONS
    elseif param.ICname != param.sourcename
        error("Source and manufactured initial conditions must match!")
    
    elseif param.ICname == "GassnerEuler"
        A = 0.1
        av = 2.0
        k = 2 * pi
        u[:,1] .= av .+ A .* sin.(k .* sum(pts, dims=2))
        @views u[:,2] .= u[:,1]
        @views u[:,3] .= u[:,1].^2

        return u
    
    elseif param.ICname == "GassnerBurgers" # BE VERY CAREFUL ABOUT MATCHING CONSTANTS
        A = 10.0
        av = 10.0
        k = 4 * pi
        u[:,1] .= A .* sin.(k .* pts) .+ av

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


    if param.sourcename == "GassnerBurgers" # ONLY IN 1D
        A = 10.0
        av = 10.0
        k = 4 * pi
        c = 2.2
        Q[:,1] .= (A * k) .* cos.(k .* (pts .- c * time)) .* (A .* sin.(k .* (pts .- c * time)) .- c .+ av)

        return Q

    elseif param.sourcename == "GassnerEuler" # CAREFUL not exactly the same as in Gassner, corrected source and solution from CPerthick (should work for 1 and 2D)
        A = 0.1
        av = 2.0
        k = 2 * pi
        c = 2.0
        Q[:,1] .=  A*k*(1.0-c) .* cos.(k * (sum(pts, dims=2) .- c*time))
        Q[:,2:end-1] .= A*k .* cos.(k .* (sum(pts, dims=2) .- c*time)) .* (2.0*av*(param.gamma-1.0)-0.5*(param.gamma-3.0)-c .+ 2.0*A*(param.gamma-1.0) .* sin.(k .* (sum(pts, dims=2) .- c*time)))
        Q[:,end] .= A*k .* cos.(k .* (sum(pts, dims=2) .- c*time)) .* (2.0*av*(param.gamma-c)-0.5*(param.gamma-1.0) .+ 2.0*A*(param.gamma-c) .* sin.(k .* (sum(pts, dims=2) .- c*time)))

        return Q
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
        else
            return make_interval(collect(range(-0.5, 0.5, param.Neldim+1)), [-1, -2])
        end
    end
end

function initialize_BCHandler(dg::DG, param::parameters)
    if param.BCname == "periodic"
        return Dict()

    elseif param.BCname == "homogeneous_Dirichlet"
        return Dict(-1 => zeros((dg.Nstates,)))
    
    elseif param.BCname == "unit_rarefaction"
        return Dict(-1 => -ones((dg.Nstates,)), -2 => ones((dg.Nstates,)))

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
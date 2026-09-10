#####################################################################
# Initial Conditions
#####################################################################

abstract type InitialCondition end

# One-state initializations

struct Sin1State <: InitialCondition
    k::Float64 # wavenumber
    phi::Float64 # phase shift
    av::Float64 # average

    function Sin1State(k::Union{Float64,Nothing}=nothing, phi::Union{Float64,Nothing}=nothing, av::Union{Float64,Nothing}=nothing)
        # These are the default values
        kdefault = 2*pi # THIS WAS CHANGED FROM 1 TO 2PI CAREFUL!
        phidefault = 0.0
        avdefault = 0.0

        new(something(k, kdefault), something(phi, phidefault), something(av, avdefault))
    end
end

function initialize_states(ic::Sin1State, dg::DG, PDE::GoverningPDE, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.Nstates != 1
        error("Sin1State only works for 1 state problems!")
    end

    u = zeros((size(pts,1), dg.Nstates))
    u .= sin.(ic.k .* (pts .- ic.phi)) .+ ic.av

    return u
end


struct Exp1State <: InitialCondition # Gaussian
    alpha::Float64 # variance for Gaussian

    function Exp1State(alpha::Union{Float64,Nothing}=nothing)
        # These are the default values
        alphadefault = -80.0

        new(something(alpha, alphadefault))
    end
end

function initialize_states(ic::Exp1State, dg::DG, PDE::GoverningPDE, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.Nstates != 1
        error("Sin1State only works for 1 state problems!")
    end

    if PDE.dim != 1 || dg.dim != 1
        error("1 State Gaussian only works in 1D!")
    end

    u = zeros((size(pts,1), dg.Nstates))
    u .= exp.(ic.alpha .* pts.^2)

    return u
end


struct Rarefaction1State <: InitialCondition # Heaviside function
    left::Float64
    right::Float64

    function Rarefaction1State(left::Union{Float64,Nothing}=nothing, right::Union{Float64,Nothing}=nothing)
        # These are the default values
        leftdefault = -1.0
        rightdefault = 1.0

        new(something(left, leftdefault), something(right, rightdefault))
    end
end

function initialize_states(ic::Rarefaction1State, dg::DG, PDE::GoverningPDE, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.Nstates != 1
        error("Sin1State only works for 1 state problems!")
    end

    if PDE.dim != 1 || dg.dim != 1
        error("1 State Rarefaction only works in 1D!")
    end

    u = ic.right .* ones(size(pts))
    u[pts .<= 0.0] .= ic.left

    return u
end


struct SuperOscillation1State <: InitialCondition # maximum oscillation at the basis nodes
    maxi::Float64
    mini::Float64

    function SuperOscillation1State(maxi::Union{Float64,Nothing}=nothing, mini::Union{Float64,Nothing}=nothing)
        # These are the default values
        maxidefault = 1.0
        minidefault = 0.0

        new(something(maxi, maxidefault), something(mini, minidefault))
    end
end

function initialize_states(ic::SuperOscillation1State, dg::DG, PDE::GoverningPDE, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.Nstates != 1
        error("SuperOscillation1State only works for 1 state problems!")
    end

    if PDE.dim != 1 || dg.dim != 1
        error("1 State Super Oscillation only works in 1D!")
    end

    u = [(ic.maxi .- ic.mini) .* Float64(isodd(i)) .+ ic.mini for i in 1:length(pts)] 
    u = reshape(u, (length(u),1))

    return u
end


struct BurgulenceIC <: InitialCondition # THIS SHOULD BE MERGED WITH BURGULENCE STRUCT AT SOME POINT
    kmax::Float64
end

function initialize_states(ic::BurgulenceIC, dg::DG, PDE::GoverningPDE, pts = nothing) # Careful, Burgulence only makes sense on unit circle!
    if !(isnothing(pts))
        error("Can only initialize Burgulence at basis nodes!")
    end

    if PDE.dim != 1 || dg.dim != 1
        error("Burgulence only works in 1D!")
    end

    Burg = Burgulence(dg.refelem.bnodestype, dg.mesh.Nel, ic.kmax)
    
    return Burg.IC
end


# Multi state initializations

struct IsentropicDensityWave <: InitialCondition # Classical test case used by (Chan 2025)
    A::Float64 # wave amplitude
    # ADD WAVE VELOCITY LATER

    function IsentropicDensityWave(A::Union{Float64,Nothing}=nothing)
        # These are the default values
        Adefault = 0.5

        new(something(A, Adefault))
    end
end

function initialize_states(ic::IsentropicDensityWave, dg::DG, PDE::EulerPerfGas, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    u = zeros((size(pts,1), dg.Nstates))
    u[:,1] .= 1 .+ ic.A .* sin.(2*pi.*sum(pts, dims=2))
    @views u[:,2] .= 0.1 .* u[:,1]

    if PDE.dim == 1
        @views u[:,end] = 10 / (PDE.gamma - 1) .+ (0.5 * 0.1^2) .* u[:,1]
    elseif PDE.dim == 2
        @views u[:,3] .= 0.2 .* u[:,1]
        @views u[:,end] = 10 / (PDE.gamma - 1) .+ (0.5 * (0.1^2 +0.2^2)) .* u[:,1]
    end

    return u
end


struct ChanWave <: InitialCondition # Scaled test case used by (Chan 2025) to investigate convergence of viscosity and entropy deficit.
    A::Float64 # wave amplitude

    function ChanWave(A::Union{Float64,Nothing}=nothing)
        # These are the default values
        Adefault = 0.5

        new(something(A, Adefault))
    end
end

function initialize_states(ic::ChanWave, dg::DG, PDE::EulerPerfGas, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.dim != 1 || dg.dim != 1
        error("Chan wave only works in 1D!")
    end

    u = zeros((size(pts,1), dg.Nstates))
    u[:,1] .= 1 .+ ic.A .* sin.(2*pi.*pts .+ 0.1/2)
    u[:,2] .= ic.A .* sin.(2*pi.*pts .+ 0.2/2) .* u[:,1]
    u[:,3] .= u[:,1].^PDE.gamma ./ (PDE.gamma-1) .+ 0.5 .* u[:,2].^2 ./ u[:,1]

    return u
end


struct GaussianVelocity <: InitialCondition # Gaussian velocity bump (ONLY IN 1D)
    rho::Float64 # constant density
    avm::Float64 # vertical bias for Gaussian
    alpham::Float64 # variance for Gaussian
    energy::Float64 # constant energy

    function GaussianVelocity(rho::Union{Float64,Nothing}=nothing, avm::Union{Float64,Nothing}=nothing, alpham::Union{Float64,Nothing}=nothing, energy::Union{Float64,Nothing}=nothing)
        # These are the default values
        rhodefault = 1.0
        avmdefault = 1.0
        alphamdefault = -80.0
        energydefault = 5.0

        new(something(rho, rhodefault), something(avm, avmdefault), something(alpham, alphamdefault), something(energy, energydefault))
    end
end

function initialize_states(ic::GaussianVelocity, dg::DG, PDE::EulerPerfGas, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.dim != 1 || dg.dim != 1
        error("Gaussian bump only works in 1D!")
    end

    u = zeros((size(pts,1), dg.Nstates))
    u[:,1] .= ic.rho
    u[:,2:end-1] .= exp.(ic.alpham .* sum(pts.^2, dims=2)) .+ ic.avm
    u[:,3] .= ic.energy

    return u
end


struct SodShockTube <: InitialCondition end

function initialize_states(ic::SodShockTube, dg::DG, PDE::EulerPerfGas, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.dim != 1 || dg.dim != 1
        error("Sod Shock Tube only works in 1D!")
    end

    u = zeros((size(pts,1), dg.Nstates))
    u[pts .< 0,1] .= 1
    u[pts .>= 0,1] .= 0.125
    u[pts .< 0,3] .= 1.0 / (PDE.gamma-1)
    u[pts .< 0,3] .= 0.1 / (PDE.gamma-1)

    return u
end


#####################################################################
# Sources and Derived initial conditions
#####################################################################

abstract type DerivedInitialCondition <: InitialCondition end

struct GassnerBurgersSource <: CLawSource
    A::Float64 # wave amplitude
    k::Float64 # wavenumber times two pi (CHANGE THIS LATER TO MAKE IT CONSISTENT WITH SIN1STATE)
    av::Float64 # average
    c::Float64 # propagation velocity

    function GassnerBurgersSource(A::Union{Float64,Nothing}=nothing, k::Union{Float64,Nothing}=nothing, av::Union{Float64,Nothing}=nothing, c::Union{Float64,Nothing}=nothing)
        # These are the default values
        Adefault = 10.0
        kdefault = 4 * pi
        avdefault = 10.0
        cdefault = 2.2

        new(something(A, Adefault), something(k, kdefault), something(av, avdefault), something(c, cdefault))
    end
end

struct GassnerBurgers <: DerivedInitialCondition
    source::GassnerBurgersSource
end

function initialize_states(ic::GassnerBurgers, dg::DG, PDE::Burgers, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.dim != 1 || dg.dim != 1
        error("GassnerBurgers only works in 1D!")
    end

    u = zeros((size(pts,1), dg.Nstates))
    u[:,1] .= ic.source.A .* sin.(ic.source.k .* pts) .+ ic.source.av

    return u
end

function compute_source(source::GassnerBurgersSource, PDE::Burgers, dg::DG, pts::Matrix{Float64}, time::Float64)
    Q = Matrix{Float64}(undef, (size(pts,1), dg.Nstates))

    if PDE.dim == 2
        println("Verifiy the accuracy of source and manufactured solution in 2D. This has not been validated...")
    end

    Q[:,1] .= (source.A * source.k) .* cos.(source.k .* (pts .- source.c * time)) .* (source.A .* sin.(source.k .* (pts .- source.c * time)) .- source.c .+ source.av)

    return Q
end


struct GassnerEulerSource <: CLawSource
    A::Float64 # wave amplitude
    k::Float64 # wavenumber times two pi (CHANGE THIS LATER TO MAKE IT CONSISTENT WITH SIN1STATE)
    av::Float64 # average
    c::Float64 # propagation velocity

    function GassnerEulerSource(A::Union{Float64,Nothing}=nothing, k::Union{Float64,Nothing}=nothing, av::Union{Float64,Nothing}=nothing, c::Union{Float64,Nothing}=nothing)
        # These are the default values
        Adefault = 0.1
        kdefault = 2 * pi
        avdefault = 2.0
        cdefault = 2.0

        new(something(A, Adefault), something(k, kdefault), something(av, avdefault), something(c, cdefault))
    end
end

struct GassnerEuler <: DerivedInitialCondition
    source::GassnerEulerSource
end

function initialize_states(ic::GassnerEuler, dg::DG, PDE::EulerPerfGas, pts = nothing)
    if isnothing(pts)
        pts = dg.bpts
    end

    if PDE.dim != 1 || dg.dim != 1
        error("GassnerEuler only works in 1D!")
    end

    u = zeros((size(pts,1), dg.Nstates))
    u[:,1] .= ic.source.av .+ ic.source.A .* sin.(ic.source.k .* sum(pts, dims=2))
    @views u[:,2] .= u[:,1]
    @views u[:,3] .= u[:,1].^2

    return u
end

function compute_source(source::GassnerEulerSource, PDE::EulerPerfGas, dg::DG, pts::Matrix{Float64}, time::Float64)
    Q = Matrix{Float64}(undef, (size(pts,1), dg.Nstates))

    if PDE.dim == 2
        println("Verifiy the accuracy of source and manufactured solution in 2D. This has not been validated...")
    end

    Q[:,1] .=  source.A*source.k*(1.0-source.c) .* cos.(source.k * (sum(pts, dims=2) .- source.c*time))
    Q[:,2:end-1] .= source.A*source.k .* cos.(source.k .* (sum(pts, dims=2) .- source.c*time)) .* (2.0*source.av*(PDE.gamma-1.0)-0.5*(PDE.gamma-3.0)-source.c .+ 2.0*source.A*(PDE.gamma-1.0) .* sin.(source.k .* (sum(pts, dims=2) .- source.c*time)))
    Q[:,end] .= source.A*source.k .* cos.(source.k .* (sum(pts, dims=2) .- source.c*time)) .* (2.0*source.av*(PDE.gamma-source.c)-0.5*(PDE.gamma-1.0) .+ 2.0*source.A*(PDE.gamma-source.c) .* sin.(source.k .* (sum(pts, dims=2) .- source.c*time)))

    return Q
end


# function initialize_states(dg::DG, param::parameters, pts = nothing)
#     flag_bpts = false

#     if isnothing(pts)
#         pts = dg.bpts
#         flag_bpts = true
#     end

#     u = zeros((size(pts,1), dg.Nstates))

#     if param.ICname == "sin_1state"
#         if param.dim != 1
#             error("IC only works in 1D!")
#         end

#         u .= sin.(2*pi .* param.k[1] .* (pts .- param.phi)) .+ param.av

#         return u
    
#     elseif param.ICname == "exp_1state"
#         if param.dim != 1
#             error("IC only works in 1D!")
#         end

#         u .= exp.(-80 .* pts.^2)

#         return u
    
#     elseif param.ICname == "rarefaction_1state"
#         if param.dim != 1
#             error("IC only works in 1D!")
#         end

#         u = ones(size(pts))
#         u[pts .<= 0.0] .= -1.0

#         return u
    
#     elseif param.ICname == "zeros_and_ones"
#         if param.dim != 1
#             error("IC only works in 1D!")
#         end

#         u = [Float64(isodd(i)) for i in 1:length(pts)] 
#         u = reshape(u, (length(u),1))

#         return u
        
#     elseif param.ICname == "IsentropicDensityWave" # Classical test case used by (Chan 2025)
#         A = 0.5 # tunable wave amplitude
#         u[:,1] .= 1 .+ A .* sin.(2*pi.*sum(pts, dims=2))
#         @views u[:,2] .= 0.1 .* u[:,1]

#         if param.dim == 1
#             @views u[:,end] = 10 / (param.gamma - 1) .+ (0.5 * 0.1^2) .* u[:,1]
#         elseif param.dim == 2
#             @views u[:,3] .= 0.2 .* u[:,1]
#             @views u[:,end] = 10 / (param.gamma - 1) .+ (0.5 * (0.1^2 +0.2^2)) .* u[:,1]
#         end

#         return u
    
#     elseif param.ICname == "ChanWave" # Scaled test case used by (Chan 2025) to investigate convergence of viscosity and entropy deficit.
#         if param.dim != 1
#             error("Only works for 1D right now!")
#         end

#         A = 0.5 # tunable wave amplitude
#         u[:,1] .= 1 .+ A .* sin.(2*pi.*pts .+ 0.1/2)
#         u[:,2] .= A .* sin.(2*pi.*pts .+ 0.2/2) .* u[:,1]
#         u[:,3] .= u[:,1].^param.gamma ./ (param.gamma-1) .+ 0.5 .* u[:,2].^2 ./ u[:,1]

#         return u
    
#     elseif param.ICname == "GaussianVelocity"
#         u[:,1] .= 1.0
#         u[:,2:end-1] .= exp.(-80 .* sum(pts.^2, dims=2)) .+ 1.0
#         u[:,3] .= 5.0

#         return u

#     elseif param.ICname == "SodShockTube"
#         if param.dim != 1
#             error("Sod shock only works in 1D!")
#         end

#         u[pts .< 0,1] .= 1
#         u[pts .>= 0,1] .= 0.125
#         u[pts .< 0,3] .= 1.0 / (param.gamma-1)
#         u[pts .< 0,3] .= 0.1 / (param.gamma-1)

#         return u
    
#     elseif param.ICname == "Burgulence"
#         if !flag_bpts
#             error("Can only initialize Burgulence at basis nodes!")
#         end

#         if (param.domain != "unit_interval_linear") || (param.BCname != "periodic")
#             error("Burgulence can only be initialized on the unit circle!")
#         end

#         bnodes = make_nodes(param.bnodes)
#         Burg = Burgulence(bnodes, param.Neldim, param.kmax)

#         return Burg.IC

    
#     # MANUFACTURED SOLUTIONS
#     elseif param.ICname != param.sourcename
#         error("Source and manufactured initial conditions must match!")
    
#     elseif param.ICname == "GassnerEuler"
#         A = 0.1
#         av = 2.0
#         k = 2 * pi
#         u[:,1] .= av .+ A .* sin.(k .* sum(pts, dims=2))
#         @views u[:,2] .= u[:,1]
#         @views u[:,3] .= u[:,1].^2

#         return u
    
#     elseif param.ICname == "GassnerBurgers" # BE VERY CAREFUL ABOUT MATCHING CONSTANTS
#         A = 10.0
#         av = 10.0
#         k = 4 * pi
#         u[:,1] .= A .* sin.(k .* pts) .+ av

#         return u
    
#     else
#         error("Unknown IC name!")
#     end
# end

# function compute_source(dg::DG, param::parameters, pts::Matrix{Float64}, time::Float64)
#     Q = Matrix{Float64}(undef, (size(pts,1), dg.Nstates))

#     if param.dim == 2
#         println("Verifiy the accuracy of source and manufactured solution in 2D. This has not been validated...")
#     end


#     if param.sourcename == "GassnerBurgers" # ONLY IN 1D
#         A = 10.0
#         av = 10.0
#         k = 4 * pi
#         c = 2.2
#         Q[:,1] .= (A * k) .* cos.(k .* (pts .- c * time)) .* (A .* sin.(k .* (pts .- c * time)) .- c .+ av)

#         return Q

#     elseif param.sourcename == "GassnerEuler" # CAREFUL not exactly the same as in Gassner, corrected source and solution from CPerthick (should work for 1 and 2D)
#         A = 0.1
#         av = 2.0
#         k = 2 * pi
#         c = 2.0
#         Q[:,1] .=  A*k*(1.0-c) .* cos.(k * (sum(pts, dims=2) .- c*time))
#         Q[:,2:end-1] .= A*k .* cos.(k .* (sum(pts, dims=2) .- c*time)) .* (2.0*av*(param.gamma-1.0)-0.5*(param.gamma-3.0)-c .+ 2.0*A*(param.gamma-1.0) .* sin.(k .* (sum(pts, dims=2) .- c*time)))
#         Q[:,end] .= A*k .* cos.(k .* (sum(pts, dims=2) .- c*time)) .* (2.0*av*(param.gamma-c)-0.5*(param.gamma-1.0) .+ 2.0*A*(param.gamma-c) .* sin.(k .* (sum(pts, dims=2) .- c*time)))

#         return Q
#     end
# end
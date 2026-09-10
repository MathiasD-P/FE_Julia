module FE_Julia
    using SparseArrays
    using LinearAlgebra
    using ForwardDiff
    using FFTW
    using FINUFFT
    using Random
    using Plots

    import FastGaussQuadrature
    
    export parameters
    export make_nodes
    export RefElemStd, RefElemSBP
    export DGStd, DGFluxDiff, DGArtVisc, DGAddRes, DGEntFilt
    export NumFlux, TPFlux, ResCorrModel, ArtViscModel, GoverningPDE
    export InitialCondition, CLawSource
    export initialize_states, initialize_mesh, initialize_BCHandler
    export build_residual!
    export ODE_solver
    export set_up_and_solve, set_up_problem

    include("mesh.jl")
    include("basis.jl")
    include("nodes.jl")
    include("refelem.jl")
    include("DG_init.jl")
    include("physics.jl")
    include("ICsandSources.jl")
    include("parameters.jl") # parameters need physics and ICs to be defined
    include("BCs.jl") # BCs need parameters to be defined
    include("build_residual.jl")
    include("ODE_solver.jl")
    include("postprocessing.jl")
    include("jacobian_residual.jl")
    include("Burgulence.jl")
    include("main.jl")
end
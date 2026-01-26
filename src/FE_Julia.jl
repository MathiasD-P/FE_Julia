module FE_Julia
    using SparseArrays
    using LinearAlgebra
    import CSV
    using Plots

    import FastGaussQuadrature
    
    export parameters
    export RefElemStd, RefElemSBP
    export DGStd, DGFluxDiff, DGArtVisc, DGAddRes, DGEntFilt
    export ODE_solver
    export set_up_and_solve

    include("mesh.jl")
    include("basis.jl")
    include("nodes.jl")
    include("refelem.jl")
    include("parameters.jl")
    include("DG_init.jl")
    include("physics.jl")
    include("BCsandICs.jl")
    include("build_residual.jl")
    include("ODE_solver.jl")
    include("postprocessing.jl")
    include("main.jl")
end
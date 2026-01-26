module FE_Julia
    export parameters
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
module FE
    export make_rectangle_quad, make_interval
    export RefElemStd
    #export empty_param
    export DGStd

    include("mesh.jl")
    include("FE_basis.jl")
    include("nodes.jl")
    include("refelem.jl")
    include("parameters.jl")
    include("DG_init.jl")
    include("physics.jl")
    include("BCsandICs.jl")
end
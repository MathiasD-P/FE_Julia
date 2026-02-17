#####################################################################
# Control parameters for the DG solver
#####################################################################

# Incomplete

mutable struct parameters
    pdetype::Union{String,Nothing}
    dim::Union{Integer,Nothing}
    dgtype::Union{String,Nothing}

    bnodes::Union{String,Nothing}
    qnodes::Union{String,Nothing}
    fnodes::Union{String,Nothing}
    enodes::Union{String,Nothing} # nodes for L2 error computation

    refelem::Union{String,Nothing}
    domain::Union{String,Nothing}
    Neldim::Union{Integer,Nothing}

    numfluxtype::Union{String,Nothing}
    twoptfluxtype::Union{String,Nothing}

    ICname::Union{String,Nothing}
    BCname::Union{String,Nothing}
    sourcename::Union{String,Nothing}
    ODE_solver::Union{String,Nothing}
    Nsteps::Union{Integer,Nothing}
    dt::Union{Real,Nothing}

    # Auxilliary computations?
    save::Bool
    calc_entropy::Bool
    OOAtest::Bool
    Nrefinements::Union{Integer,Nothing}

    # And all other physical constants for the problem...
    gamma # specific heat ratio for Euler
    a # advection speed for lin advection

    function parameters(;
                     pdetype=nothing,
                     dim=nothing,
                     dgtype=nothing,
                     bnodes=nothing,
                     qnodes=nothing,
                     enodes=nothing,
                     fnodes=nothing,
                     refelem=nothing,
                     domain=nothing,
                     Neldim=nothing,
                     numfluxtype=nothing,
                     twoptfluxtype=nothing,
                     ICname=nothing,
                     BCname=nothing,
                     sourcename=nothing,
                     ODE_solver=nothing,
                     Nsteps=nothing,
                     dt=nothing,
                     save=false,
                     calc_entropy=false,
                     OOAtest=false,
                     Nrefinements=nothing,
                     gamma=nothing,
                     a=nothing)

                     new(pdetype, dim, dgtype, bnodes, qnodes, fnodes, enodes, refelem, domain, Neldim, numfluxtype, twoptfluxtype, ICname, BCname, sourcename, ODE_solver, Nsteps, dt, save, calc_entropy, OOAtest, Nrefinements, gamma,a)
    end
end

function parse_parameters(param_file)
    return
end


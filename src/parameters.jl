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
    qmnodes::Union{String,Nothing}
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

    AVcoeff::Union{String,Nothing} # only for artificial viscosity
    Rescorr::Union{String,Nothing} # only for AddRes

    # Auxilliary computations?
    save::Bool
    calc_entropy::Bool
    OOAtest::Bool
    Nrefinements::Union{Integer,Nothing}
    dtlim::Union{Vector{Real}, Nothing}
    maxval::Union{Real,Nothing}

    # And all other physical constants for the problem...
    gamma # specific heat ratio for Euler
    a # advection speed for lin advection
    addviscosity # viscosity offset for artificial viscosity
    k::Union{Vector{Real},Nothing} # wavenumber for sinusoidal initializations (by state)
    phi::Union{Vector{Real},Nothing} # phase shift for sinusoidal initializations (by state)
    av::Union{Vector{Real},Nothing} # average for sinusoidal initializations (by state)
    kmax # wavenumber cutoff for Burgulence

    function parameters(;
                     pdetype=nothing,
                     dim=nothing,
                     dgtype=nothing,

                     bnodes=nothing,
                     qnodes=nothing,
                     qmnodes=nothing,
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

                     AVcoeff=nothing,
                     Rescorr=nothing,

                     save=false,
                     calc_entropy=false,
                     OOAtest=false,
                     Nrefinements=nothing,
                     dtlim=nothing,
                     maxval=nothing,

                     gamma=nothing,
                     a=nothing,
                     addviscosity=nothing,
                     k=nothing,
                     phi=nothing,
                     av=nothing,
                     kmax=nothing)

                     new(pdetype,
                         dim,
                         dgtype,
                         bnodes,
                         qnodes,
                         qmnodes,
                         fnodes,
                         enodes,
                         refelem,
                         domain,
                         Neldim,
                         numfluxtype,
                         twoptfluxtype,
                         ICname,
                         BCname,
                         sourcename,
                         ODE_solver,
                         Nsteps,
                         dt,
                         AVcoeff,
                         Rescorr,
                         save,
                         calc_entropy,
                         OOAtest,
                         Nrefinements,
                         dtlim,
                         maxval,
                         gamma,
                         a,
                         addviscosity,
                         k,
                         phi,
                         av,
                         kmax)
    end
end

function parse_parameters(param_file)
    return
end


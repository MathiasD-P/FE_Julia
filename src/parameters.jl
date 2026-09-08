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

function parse_parameters_PDE(param::parameters)
    if param.pdetype == "LinAdv"
        PDE = PDE_LinAdv(param.dim, param.a)
    elseif param.pdetype == "Burgers"
        PDE = PDE_Burgers(param.dim)
    elseif param.pdetype == "EulerPerfGas"
        PDE = PDE_EulerPerfGas(param.dim, param.gamma)
    else
        error("The PDE type $(param.pdetype) is not defined. Please choose from LinAdv, Burgers, or EulerPerfGas.")
    end

    return PDE
end

function parse_parameters_numflux(param::parameters)
    if param.numfluxtype == "central"
        numflux = CentralNumFlux()
    elseif param.numfluxtype == "upwind"
        numflux = UpwindNumFlux()
    elseif param.numfluxtype == "EC_split"
        numflux = ECSplitNumFlux()
    elseif param.numfluxtype == "LF"
        numflux = LFNumFlux()
    elseif param.numfluxtype == "EC_Chandrashekar"
        numflux = ECChandrashekarNumFlux()
    elseif param.numfluxtype == "ES_Chandrashekar_dissip"
        numflux = ESChandrashekarDissipNumFlux()
    else
        error("The numerical flux $(param.numfluxtype) is not defined.")
    end

    return numflux
end

function parse_parameters_tpflux(param::parameters)
    if param.twoptfluxtype == "EC_split"
        tpflux = ECSplitTPFlux()
    elseif param.twoptfluxtype == "AV_split"
        tpflux = AVSplitTPFlux()
    elseif param.twoptfluxtype == "OGL_split"
        if occursin(r"^\(.*\)-.*$", param.qnodes) # Check if this is the name of a tensor prod node
            m = match(r"\((.*)\)-(.*)", param.qnodes)
            q = parse.(Int, split(m.captures[1], 'x')) - 1
        else
            error("Cannot create an OGL split two poin flux from the specified quadrature node name!")
        end
        alpha = (q+3) / (2*q+3)
        beta = q / (2*q+3)
        tpflux = SplitTPFlux(alpha, beta)
    elseif param.twoptfluxtype == "OGLL_split"
        if occursin(r"^\(.*\)-.*$", param.qnodes) # Check if this is the name of a tensor prod node
            m = match(r"\((.*)\)-(.*)", param.qnodes)
            q = parse.(Int, split(m.captures[1], 'x')) - 1
        else
            error("Cannot create an OGLL split two poin flux from the specified quadrature node name!")
        end
        alpha = (q+1) / (2*q+1)
        beta = q / (2*q+1)
        tpflux = SplitTPFlux(alpha, beta)
    else
        error("The two point flux $(param.twoptfluxtype) is not defined.")
    end

    return tpflux
end

function parse_parameters_artvisc(param::parameters)
    if param.dgtype != "DGArtVisc"
        println("Warning: You are using the artificial viscosity parameter parser for a non-compatible DG type!")
    end

    if param.AVcoeff == "AVdissip"
        artviscmodel = AVdissip(param.addviscosity)
    elseif param.AVcoeff == "AVEC"
        artviscmodel = AVEC(param.addviscosity)
    elseif param.AVcoeff == "NoAV"
        artviscmodel = NoAV(param.addviscosity)
    else
        error("The artificial viscosity model $(param.AVcoeff) is not defined.")
    end

    return artviscmodel
end

function parse_parameters_rescorr(param::parameters)
    if param.dgtype != "DGAddRes"
        println("Warning: You are using the residual correction parameter parser for a non-compatible DG type!")
    end

    if param.Rescorr == "Rescorrdissip"
        rescorrmodel = ResCorrDissip()
    elseif param.Rescorr == "RescorrEC"
        rescorrmodel = ResCorrEC()
    elseif param.Rescorr == "NoRescorr"
        rescorrmodel = NoResCorr()
    else
        error("The residual correction model $(param.Rescorr) is not defined.")
    end

    return rescorrmodel
end
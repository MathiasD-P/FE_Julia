#####################################################################
# Evaluate the jacobian of the residual for our DG implementations
#####################################################################

function residual_jacobian(u, t::Float64, dg::DGStd, physics::PhysProp)
    # First, we wrap the build_residual! function
    function residual_wrapper(u)
        residual = similar(u)
        build_residual!(residual, u, t, dg, physics)
        return residual
    end

    return ForwardDiff.jacobian(residual_wrapper, u)
end

function residual_jacobian_Dirichlet(u, t::Float64, dg::Union{DGStd, DGFluxDiff, DGAddRes}, physics::PhysProp)
    if dg.dim !=1 || physics.PDE.dim != 1
        error("residual_jacobian_Dirichlet only works in one dimension.")
    end

    function residual_wrapper(uext)
        newphysics = deepcopy(physics)
        newphysics.BChandler = Dict(-1 => uext[1,:], -2 => uext[2,:])
        residual = similar(uext, size(uext,1)-2, physics.PDE.Nstates)
        build_residual!(residual, uext[3:end,:], t, dg, newphysics)

        return residual
    end

    uext = [reshape(physics.BChandler[-1], 1, physics.PDE.Nstates) ; reshape(physics.BChandler[-2], 1, physics.PDE.Nstates) ; u]

    return ForwardDiff.jacobian(residual_wrapper, uext)
end
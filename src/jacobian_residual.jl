#####################################################################
# Evaluate the jacobian of the residual for our DG implementations
#####################################################################

function residual_jacobian(u, t::Float64, BChandler::Dict, dg::DGStd, param::parameters)
    # First, we wrap the build_residual! function
    function residual_wrapper(u)
        residual = similar(u)
        build_residual!(residual, u, t, BChandler, dg, param)
        return residual
    end

    return ForwardDiff.jacobian(residual_wrapper, u)
end
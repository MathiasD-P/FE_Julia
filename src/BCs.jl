#####################################################################
# Consistent mesh and BCHandler initializations (must absolutely be updated together)
#####################################################################

function initialize_mesh(param::parameters)
    if param.domain == "unit_square_linear_quad"
        if param.BCname == "periodic"
            return make_rectangle_quad(param.Neldim, param.Neldim)
        elseif param.BCname == "homogeneous_Dirichlet"
            return make_rectangle_quad(param.Neldim, param.Neldim, [-1,-1,-1,-1])
        end

    elseif param.domain == "unit_interval_linear"
        if param.BCname == "periodic"
            return make_interval(collect(range(-0.5, 0.5, param.Neldim+1)), [0, 0])
        elseif param.BCname == "homogeneous_Dirichlet"
            return make_interval(collect(range(-0.5, 0.5, param.Neldim+1)), [-1, -1])
        else
            return make_interval(collect(range(-0.5, 0.5, param.Neldim+1)), [-1, -2])
        end
    end
end

function initialize_BCHandler(dg::DG, param::parameters)
    if param.BCname == "periodic"
        return Dict()

    elseif param.BCname == "homogeneous_Dirichlet"
        return Dict(-1 => zeros((dg.Nstates,)))
    
    elseif param.BCname == "unit_rarefaction"
        return Dict(-1 => -ones((dg.Nstates,)), -2 => ones((dg.Nstates,)))

    elseif param.BCname == "SodShockTube"
        return Dict(-1 => [1.0, 0.0, 1.0 / (param.gamma-1)], -2 => [0.125, 0.0, 0.1 / (param.gamma-1)])
    end
end

function evaluate_BC(BCHandler::Dict, dg::DG, t)
    if isempty(BCHandler)
        return spzeros(dg.NFval, dg.Nstates)

    elseif BCHandler isa Dict{<:Integer, <:Vector{<:Real}}
        return sum([dg.BFtoF[itag] * BCHandler[dg.mesh.BCtags[itag]]' for itag in 1:dg.mesh.Ntags])
    end
end

function isperiodic(BCHandler::Dict{Integer, Any})
    return isempty(BCHandler)
end

function isperiodicunitinterval(mesh::LMesh)
    check_periodic = sum(mesh.BCtags) == 0
    check_dim = mesh.dim == 1
    check_length = abs(maximum(maximum, mesh.elements) - minimum(minimum, mesh.elements) - 1.0) < 1e-12

    return check_dim && check_periodic && check_length
end
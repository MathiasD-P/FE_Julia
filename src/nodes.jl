#####################################################################
# Create, manipulate and store volume and face nodes
# -> Very dirty code, need to rewrite extract_face_operators and extract_volume_operators as a single function
# -> Not really optimized, but shouldn't be a too big bottleneck since most arrays are evaluated in
#    the init phase.
# AT SOME POINT, WILL HAVE TO SWITCH TO LAZY ALLOCATIONS OF
#####################################################################

abstract type AbstractNodes end

mutable struct Tensorprod_nodes <: AbstractNodes
    name::String # name of the nodes
    K::Integer # number of nodes in each direction
    dim::Integer # dimension of space
end

mutable struct Fnodes
    refshape::String # specifies N immersions from the node space to the ambient space
    nodes::AbstractNodes # nodes of dim-1
end

function make_nodes(name::String) # Create nodes from their name
    if occursin(r"^\(.*\)-.*$", name) # Check if this is the name of a tensor prod node
        m = match(r"\((.*)\)-(.*)", name)
        K = parse.(Int, split(m.captures[1], 'x'))
        dim = count(==('x'), m.captures[1]) + 1

        if all(K.==K[1])
            K = K[1]
        else
            error("Number of nodes in each spatial direction must agree for tensor prod nodes!")
        end

        return Tensorprod_nodes(m.captures[2], K, dim)
    end

    error("Unknown node name!")
end

function make_nodes(refshape::String, name::String) # Create Fnodes from their name
    return Fnodes(refshape, make_nodes(name))
end

function evaluate(nodes::Tensorprod_nodes)
    if nodes.name == "GLL"
        xi, _ = FastGaussQuadrature.gausslobatto(nodes.K)
    elseif nodes.name == "GL"
        xi, _ = FastGaussQuadrature.gaussjacobi(nodes.K,0.0,0.0)
    elseif nodes.name == "linear"
        xi = collect(range(0, 1, K))
    end
    xi = 0.5 * (reshape(xi, length(xi), 1) .+ 1)

    if nodes.dim == 0
        return reshape([0],1,1)
    elseif nodes.dim == 1
        return xi
    elseif nodes.dim == 2
        return point_tensorprod(xi,xi)
    end
end

function evaluate(fnodes::Fnodes)
    if fnodes.refshape == "interval"
        return [reshape([0],1,1), reshape([1],1,1)]
    elseif fnodes.refshape == "square"
        xi = evaluate(fnodes.nodes)
        return [[xi zeros(size(xi))], [ones(size(xi)) xi], [reverse(xi) ones(size(xi))], [zeros(size(xi)) reverse(xi)]]
    end
end

function numnodes(nodes::Tensorprod_nodes)
    return nodes.K^nodes.dim
end

function numfaces(fnodes::Fnodes)
    if fnodes.refshape == "square"
        return 4
    elseif fnodes.refshape == "interval"
        return 2
    end
end

function Base.show(io::IO, nodes::Tensorprod_nodes)
    message = "("
    for i = 1:nodes.dim
        message *= string(nodes.K)
        if i != nodes.dim
            message *= "x"
        end
    end
    message *= ")-" * nodes.name
    print(io,message)
end

function extract_volume_quadrature(basisnodes::Tensorprod_nodes, volnodes::Tensorprod_nodes)
    if basisnodes.dim != volnodes.dim
        throw(ArgumentError("Dimension of basis and volume nodes must agree!"))
    end

    if volnodes.name == "GLL"
        xi_vol, w = FastGaussQuadrature.gausslobatto(volnodes.K)
    elseif volnodes.name == "GL"
        xi_vol, w = FastGaussQuadrature.gaussjacobi(volnodes.K,0.0,0.0)
    end
    xi_vol = 0.5 * (reshape(xi_vol, length(xi_vol), 1) .+ 1)
    w *= 0.5

    if (basisnodes.name == volnodes.name) & (basisnodes.K == volnodes.K)
        xi_basis = xi_vol
        V = Matrix(I(basisnodes.K)).+0.0
    else
        if basisnodes.name == "GLL"
            xi_basis, _ = FastGaussQuadrature.gausslobatto(basisnodes.K)
        elseif basisnodes.name == "GL"
            xi_basis, _ = FastGaussQuadrature.gaussjacobi(basisnodes.K,0.0,0.0)
        end
        xi_basis = 0.5 * (reshape(xi_basis, length(xi_basis), 1) .+ 1)

        V = Lagrange_Vandermonde1D(xi_vol, xi_basis)
    end

    if basisnodes.dim == 1
        return V, vec(w), xi_vol
    elseif basisnodes.dim == 2
        chivol = Vandermonde_tensorprod(V, V)
        wvol = vec(w * transpose(w))
        pts = point_tensorprod(xi_vol,xi_vol)
        return chivol, wvol, pts
    end
end

function extract_volume_operators(basisnodes::Tensorprod_nodes, volnodes::Tensorprod_nodes)
    if basisnodes.dim != volnodes.dim
        throw(ArgumentError("Dimension of basis and volume nodes must agree!"))
    end

    if volnodes.name == "GLL"
        xi_vol, w = FastGaussQuadrature.gausslobatto(volnodes.K)
    elseif volnodes.name == "GL"
        xi_vol, w = FastGaussQuadrature.gaussjacobi(volnodes.K,0.0,0.0)
    end
    xi_vol = 0.5 * (reshape(xi_vol, length(xi_vol), 1) .+ 1)
    w *= 0.5

    if (basisnodes.name == volnodes.name) & (basisnodes.K == volnodes.K)
        xi_basis = xi_vol
        V = Matrix(I(basisnodes.K)).+0.0
    else
        if basisnodes.name == "GLL"
            xi_basis, _ = FastGaussQuadrature.gausslobatto(basisnodes.K)
        elseif basisnodes.name == "GL"
            xi_basis, _ = FastGaussQuadrature.gaussjacobi(basisnodes.K,0.0,0.0)
        end
        xi_basis = 0.5 * (reshape(xi_basis, length(xi_basis), 1) .+ 1)

        V = Lagrange_Vandermonde1D(xi_vol, xi_basis)
    end

    dV = DLagrange_Vandermonde1D(xi_vol, xi_basis)

    if basisnodes.dim == 1
        return V, vec(w), [dV], xi_basis
    elseif basisnodes.dim == 2
        chivol = Vandermonde_tensorprod(V, V)
        wvol = vec(w * transpose(w))
        dchivol = [Vandermonde_tensorprod(dV, V), Vandermonde_tensorprod(V, dV)]
        solpts = point_tensorprod(xi_basis,xi_basis)
        return chivol, wvol, dchivol, solpts
    end
end

function extract_face_operators(basisnodes::AbstractNodes, facenodes::Fnodes)
    if basisnodes.name == "GLL"
            xi_basis, _ = FastGaussQuadrature.gausslobatto(basisnodes.K)
        elseif basisnodes.name == "GL"
            xi_basis, _ = FastGaussQuadrature.gaussjacobi(basisnodes.K,0.0,0.0)
        end
        xi_basis = 0.5 * (reshape(xi_basis, length(xi_basis), 1) .+ 1)

    if facenodes.refshape == "interval"
        chiface = [Lagrange_Vandermonde1D([0], xi_basis), Lagrange_Vandermonde1D([1], xi_basis)]
        facepts = [reshape([0],1,1), reshape([1],1,1)]
        w = [[1], [1]]

        return chiface, w, facepts

    elseif facenodes.refshape == "square"
        if facenodes.nodes.name == "GLL"
            xi_face, w = FastGaussQuadrature.gausslobatto(facenodes.nodes.K)
        elseif facenodes.nodes.name == "GL"
            xi_face, w = FastGaussQuadrature.gaussjacobi(facenodes.nodes.K,0.0,0.0)
        end

        xi_face = 0.5 * (reshape(xi_face, length(xi_face), 1) .+ 1)
        w = 0.5 * [vec(w), vec(w), reverse(vec(w)), reverse(vec(w))]
        facepts = [[xi_face zeros(size(xi_face))], [ones(size(xi_face)) xi_face], [reverse(xi_face) ones(size(xi_face))], [zeros(size(xi_face)) reverse(xi_face)]]

        Vxi = Lagrange_Vandermonde1D(xi_face, xi_basis)
        V0 = Lagrange_Vandermonde1D([0], xi_basis)
        V1 = Lagrange_Vandermonde1D([1], xi_basis)

        chiface = [Vandermonde_tensorprod(Vxi, V0),
                   Vandermonde_tensorprod(V1, Vxi),
                   Vandermonde_tensorprod(reverse(Vxi,dims=1), V1),
                   Vandermonde_tensorprod(V0, reverse(Vxi, dims=1))]
        return chiface, w, facepts
    end
end

function select_exterior_mtrx(mesh::AbstractMesh, Nfnodes::Integer)
    face_mtrx = spzeros(mesh.Nel * mesh.Nface * Nfnodes, mesh.Nel * mesh.Nface * Nfnodes)
    for ielem in 1:mesh.Nel
        for iface in 1:mesh.Nface
            if mesh.connectivity[ielem,iface][1] > 0 # BOUNDARIES ARE IGNORED!
                i1 = (ielem-1) * mesh.Nface * Nfnodes + (iface-1) * Nfnodes + 1
                i2 = (ielem-1) * mesh.Nface * Nfnodes + iface * Nfnodes
                j1 = (mesh.connectivity[ielem,iface][1]-1) * mesh.Nface * Nfnodes + (mesh.connectivity[ielem,iface][2]-1) * Nfnodes + 1
                j2 = (mesh.connectivity[ielem,iface][1]-1) * mesh.Nface * Nfnodes + mesh.connectivity[ielem,iface][2] * Nfnodes
                face_mtrx[CartesianIndex.(i1:i2, j2:-1:j1)] = ones(Nfnodes)
            end
        end
    end
    return face_mtrx
end

function select_BC_mtrx(mesh::AbstractMesh, Nfnodes::Integer, BC::String)
    BC_mtrx = spzeros(mesh.Nel * mesh.Nface * Nfnodes, mesh.Nel * mesh.Nface * Nfnodes)
    for ielem in 1:mesh.Nel
        for iface in 1:mesh.Nface
            if (mesh.connectivity[ielem,iface][2] == 0) & (BC == "Dirichlet")
                i1 = (ielem-1) * mesh.Nface * Nfnodes + (iface-1) * Nfnodes + 1
                i2 = (ielem-1) * mesh.Nface * Nfnodes + iface * Nfnodes
                BC_mtrx[CartesianIndex.(i1:i2, i1:i2)] = ones(Nfnodes)
            elseif (mesh.connectivity[ielem,iface][2] == -1) & (BC == "Neumann")
                i1 = (ielem-1) * mesh.Nface * Nfnodes + (iface-1) * Nfnodes + 1
                i2 = (ielem-1) * mesh.Nface * Nfnodes + iface * Nfnodes
                BC_mtrx[CartesianIndex.(i1:i2, i1:i2)] = ones(Nfnodes)
            end
        end
    end
    return BC_mtrx
end
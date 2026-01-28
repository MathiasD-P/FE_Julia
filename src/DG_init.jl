abstract type DG end

mutable struct DGStd <: DG
    DOF::Int64
    NFval::Int64
    Nstates::Int64
    dim::Int64

    mesh::AbstractMesh
    refelem::RefElemStd
    nphys::Matrix{Float64} # unit physical normals at all the face points

    FtoF::SparseMatrixCSC{Float64, Int64} # (NFval, NFval) Maps globally ordered face values to globally ordered connected face values
    BFtoF::Vector{SparseVector{Float64, Int64}} # (Ntags,) Vector  of (,NFval) vectors identifying the face nodes for each tag

    bpts::Matrix{Float64} # (DOF, dim) matrix storing the physical coordinates of the basis points
    bfpts::Vector{Matrix{Float64}} # (Ntags,) Vector of (Ntaggedpts, dim) vectors for the physical coordinates of the face points tagged

    function DGStd(Nstates::Integer, refelem::RefElemStd, mesh::AbstractMesh)
        DOF = refelem.Nbnodes * mesh.Nel
        NFval = refelem.Nfnodes * refelem.Nfaces * mesh.Nel

        if mesh.dim != refelem.dim
            error("dimension of mesh and refelem must agree!")
        else
            dim = mesh.dim
        end

        nphys = compute_unitphys_normals(mesh, refelem)

        FtoF = compute_FtoF(mesh, refelem)
        BFtoF = compute_BFtoF(mesh, refelem)

        bpts = reduce(vcat, mapping(mesh, refelem.bnodes, ielem) for ielem in 1:mesh.Nel)
        fpts = reduce(vcat, mapping(mesh, refelem.fnodes, ielem) for ielem in 1:mesh.Nel)
        bfpts = [fpts[BFtoF[itag] .== 1,:] for itag in 1:mesh.Ntags]

        new(
            DOF,
            NFval,
            Nstates,
            dim,

            mesh,
            refelem,
            nphys,

            FtoF,
            BFtoF,

            bpts,
            bfpts
        )
    end
end

mutable struct DGFluxDiff <: DG
    DOF::Int64
    NFval::Int64
    Nstates::Int64
    dim::Int64

    mesh::AbstractMesh
    refelem::RefElemSBP
    nphys::Matrix{Float64} # unit physical normals at all the face points

    FtoF::SparseMatrixCSC{Float64, Int64} # (NFval, NFval) Maps globally ordered face values to globally ordered connected face values
    BFtoF::Vector{SparseVector{Float64, Int64}} # (Ntags,) Vector  of (,NFval) vectors identifying the face nodes for each tag

    bpts::Matrix{Float64} # (DOF, dim) matrix storing the physical coordinates of the basis points
    bfpts::Vector{Matrix{Float64}} # (Ntags,) Vector of (Ntaggedpts, dim) vectors for the physical coordinates of the face points tagged

    function DGFluxDiff(Nstates::Integer, refelem::RefElemSBP, mesh::AbstractMesh)
        DOF = refelem.Nbnodes * mesh.Nel
        NFval = refelem.Nfnodes * refelem.Nfaces * mesh.Nel

        if mesh.dim != refelem.dim
            error("dimension of mesh and refelem must agree!")
        else
            dim = mesh.dim
        end

        nphys = compute_unitphys_normals(mesh, refelem)

        FtoF = compute_FtoF(mesh, refelem)
        BFtoF = compute_BFtoF(mesh, refelem)

        bpts = reduce(vcat, mapping(mesh, refelem.bnodes, ielem) for ielem in 1:mesh.Nel)
        fpts = reduce(vcat, mapping(mesh, refelem.fnodes, ielem) for ielem in 1:mesh.Nel)
        bfpts = [fpts[BFtoF[itag] .== 1,:] for itag in 1:mesh.Ntags]

        new(
            DOF,
            NFval,
            Nstates,
            dim,

            mesh,
            refelem,
            nphys,

            FtoF,
            BFtoF,

            bpts,
            bfpts
        )
    end
end

mutable struct DGArtVisc <: DG
    DOF::Integer
end

mutable struct DGAddRes <: DG
    DOF::Integer
end

mutable struct DGEntFilt <: DG
    DOF::Integer
end

function compute_unitphys_normals(mesh::LMesh, refelem::AbstractRefElem)
    if all(x -> x === mesh.CT[1], mesh.CT) # check if cofactor matrix is constant over all elements
        unitn = vcat((n -> mesh.CT[1]' * n ./ norm(n)).(refelem.nref)...)

        compute = false
    else
        compute = true
    end

    allnormals = Matrix{Float64}(undef, refelem.Nfnodes * refelem.Nfaces * mesh.Nel, mesh.dim)

    for ielem in 1:mesh.Nel
        if compute
            unitn = vcat((n -> mesh.CT[ielem]' * n ./ norm(n)).(refelem.nref)...)
        end

        for iface in 1:refelem.Nfaces
            index = (ielem-1)*refelem.Nfnodes*refelem.Nfaces + (iface-1)*refelem.Nfnodes + 1
            allnormals[index:index+refelem.Nfnodes-1,:] .= unitn[iface,:]
        end
    end

    return allnormals
end

function compute_FtoF(mesh::AbstractMesh, refelem::AbstractRefElem)
    face_mtrx = spzeros(mesh.Nel * refelem.Nfaces * refelem.Nfnodes, mesh.Nel * refelem.Nfaces * refelem.Nfnodes)

    for ielem in 1:mesh.Nel
        for iface in 1:refelem.Nfaces
            if mesh.connectivity[ielem,iface][1] > 0 # BOUNDARIES ARE IGNORED!
                i1 = (ielem-1) * refelem.Nfaces * refelem.Nfnodes + (iface-1) * refelem.Nfnodes + 1
                i2 = (ielem-1) * refelem.Nfaces * refelem.Nfnodes + iface * refelem.Nfnodes
                j1 = (mesh.connectivity[ielem,iface][1]-1) * refelem.Nfaces * refelem.Nfnodes + (mesh.connectivity[ielem,iface][2]-1) * refelem.Nfnodes + 1
                j2 = (mesh.connectivity[ielem,iface][1]-1) * refelem.Nfaces * refelem.Nfnodes + mesh.connectivity[ielem,iface][2] * refelem.Nfnodes
                face_mtrx[CartesianIndex.(i1:i2, j2:-1:j1)] = ones(refelem.Nfnodes)
            end
        end
    end

    return face_mtrx
end

function compute_BFtoF(mesh::AbstractMesh, refelem::AbstractRefElem)
    function compute_BFtoF_tag(mesh::AbstractMesh, refelem::AbstractRefElem, BCtag::Integer)
        BC_mtrx = spzeros(mesh.Nel * refelem.Nfaces * refelem.Nfnodes,)
        for ielem in 1:mesh.Nel
            for iface in 1:mesh.Nface
                if (mesh.connectivity[ielem,iface][2] == BCtag)
                    i1 = (ielem-1) * mesh.Nface * refelem.Nfnodes + (iface-1) * refelem.Nfnodes + 1
                    i2 = (ielem-1) * mesh.Nface * refelem.Nfnodes + iface * refelem.Nfnodes
                    BC_mtrx[i1:i2] = ones(refelem.Nfnodes)
                end
            end
        end
        return BC_mtrx
    end

    return [compute_BFtoF_tag(mesh, refelem, mesh.BCtags[itag]) for itag in 1:mesh.Ntags]
end
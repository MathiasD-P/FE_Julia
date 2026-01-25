# NOTE: CHANGE PACKAGING OF SURFACE TERMS EVENTUALLY.

using SparseArrays

abstract type AbstractRefElem end

struct RefElemStd <: AbstractRefElem
    shape::String # type of element

    dim::Integer # dimension of element
    Nfaces::Integer
    Nbnodes::Integer # total number of basis nodes
    Nqnodes::Integer # total number of quadrature nodes
    Nfnodes::Integer # number of face nodes on EACH face

    bnodestype::AbstractNodes
    qnodestype::AbstractNodes
    fnodestype::Fnodes

    bnodes::Matrix{Float64}
    fnodes::Matrix{Float64}

    nref::AbstractArray # reference normal vectors

    chiq::AbstractArray # volume quadrature Vandermonde matrix
    chif::AbstractArray # face quadrature Vandermonde matrix
    wq::AbstractArray # volume quadrature weights
    wf::AbstractArray # face quadrature weights

    M::AbstractArray # reference mass matrix
    bh::AbstractArray # boundary integration weights times ref normal
    Dh::AbstractArray # reference differentiation matrix in each direction
    Ph::AbstractArray # reference projection matrix
    LIFT::AbstractArray # lifting operator in each direction

    function RefElemStd(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes)
        shape = fnodestype.refshape

        dim = bnodestype.dim
        Nfaces = numfaces(fnodestype)
        Nbnodes = numnodes(bnodestype)
        Nqnodes = numnodes(qnodestype)
        Nfnodes = numnodes(fnodestype.nodes)

        nref = unitnormals(fnodestype)

        bnodestype = bnodestype
        qnodestype = qnodestype
        fnodestype = fnodestype

        chiq, wq, dchi, bnodes = extract_volume_operators(bnodestype, qnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        M = chiq' * Diagonal(wq) * chiq
        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([sparse(chif[iface]) for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]
        Ph = sparse(M \ chiq' * Diagonal(wq))
        Dh = (d -> M \ (chiq' * Diagonal(wq)) * d).(dchi)
        LIFT = (b -> (M \ chif') * Diagonal(b)).(bh)

        new(
            shape,

            dim,
            Nfaces,
            Nbnodes,
            Nqnodes,
            Nfnodes,

            bnodestype,
            qnodestype,
            fnodestype,

            bnodes,
            fnodes,

            nref,

            chiq,
            chif,
            wq,
            wf,

            M,
            bh,
            Dh,
            Ph,
            LIFT
        )
    end
end

function Base.show(io::IO, re::AbstractRefElem)
    print(io, "Reference element type: " * re.shape * "\n")
    print(io, "Basis nodes: \n")
    show(io, re.bnodestype)
    print(io,"\nVolume quadrature nodes: \n")
    show(io, re.qnodestype)
    print(io,"\nFace quadrature nodes: \n")
    show(io, re.fnodestype.nodes)
    print(io,"\nI look like this: \n\n")
    
    if re.shape == "square"
        print(io,"(0,1) --[f3]-- (1,1) \n  |  (n2)  (n4)  | \n [f4]           [f2]\n  |  (n1)  (n3)  | \n(0,0) --[f1]-- (1,0)\n")
    elseif re.shape == "interval"
        print(io,"(0)[f1]---(n1)---(n2)---[f2](1)\n")
    end
end

function unitnormals(fnodes::Fnodes)
    if fnodes.refshape == "square"
        return [[0.0, -1.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]
    elseif fnodes.refshape == "interval"
        return [[-1.0], [1.0]]
    end
end
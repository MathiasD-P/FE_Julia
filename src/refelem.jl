# NOTE: CHANGE PACKAGING OF SURFACE TERMS EVENTUALLY.

abstract type AbstractRefElem end

struct RefElemStd <: AbstractRefElem
    shape::String # type of element

    dim::Integer # dimension of element
    Nfaces::Integer # total number of faces
    Nbnodes::Integer # total number of basis nodes
    Nqnodes::Integer # total number of quadrature nodes
    Nfnodes::Integer # number of face nodes on EACH face

    bnodestype::AbstractNodes
    qnodestype::AbstractNodes
    fnodestype::Fnodes

    bnodes::Matrix{Float64}
    qnodes::Matrix{Float64}
    fnodes::Matrix{Float64}

    nref::Vector{Vector{Float64}} # reference normal vectors

    chiq::Matrix{Float64} # volume quadrature Vandermonde matrix
    chif::Matrix{Float64} # face quadrature Vandermonde matrix
    wq::Vector{Float64} # volume quadrature weights

    M::Matrix{Float64}# reference mass matrix
    bh::Vector{Vector{Float64}} # boundary integration weights times ref normal
    Dh::Vector{Matrix{Float64}} # reference differentiation matrix in each direction
    Qh::Vector{Matrix{Float64}} # reference stiffness matrix in each direction
    MinvQhT::Vector{Matrix{Float64}} # for weak DG on linear meshes
    Ph::Matrix{Float64} # reference projection matrix
    LIFT::Vector{Matrix{Float64}} # lifting operator in each direction

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

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Ph = (M \ chiq' * Diagonal(wq))
        Qh = (d -> chiq' * Diagonal(wq) * d).(dchi)
        Dh = (q -> M \ q).(Qh)
        MinvQhT = (q -> M \ q').(Qh)
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
            qnodes,
            fnodes,

            nref,

            chiq,
            chif,
            wq,

            M,
            bh,
            Dh,
            Qh,
            MinvQhT,
            Ph,
            LIFT
        )
    end

    function RefElemStd(bnodestype::AbstractNodes, qnodestype::AbstractNodes, qmnodestype::AbtractNodes, fnodestype::Fnodes)
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

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chimq, wmq, _, _, _ = extract_volume_operators(bnodestype, qmnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Mm = chimq' * Diagonal(wmq) * chimq
        Ph = (M \ chiq' * Diagonal(wq))
        Qh = (d -> chiq' * Diagonal(wq) * d).(dchi)
        Dh = (q -> Mm \ q).(Qh)
        MinvQhT = (q -> M \ q').(Qh)
        LIFT = (b -> (Mm \ chif') * Diagonal(b)).(bh)

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
            qnodes,
            fnodes,

            nref,

            chiq,
            chif,
            wq,

            M,
            bh,
            Dh,
            Qh,
            MinvQhT,
            Ph,
            LIFT
        )
    end
end

struct RefElemSBP <: AbstractRefElem
    shape::String # type of element

    dim::Integer # dimension of element
    Nfaces::Integer # total number of faces
    Nbnodes::Integer # total number of basis nodes
    Nqnodes::Integer # total number of quadrature nodes
    Nfnodes::Integer # number of face nodes on EACH face

    bnodestype::AbstractNodes
    qnodestype::AbstractNodes
    fnodestype::Fnodes

    bnodes::Matrix{Float64}
    qnodes::Matrix{Float64}
    fnodes::Matrix{Float64}

    nref::Vector{Vector{Float64}} # reference normal vectors

    chiq::Matrix{Float64} # volume quadrature Vandermonde matrix
    chif::Matrix{Float64} # face quadrature Vandermonde matrix
    wq::Vector{Float64} # volume quadrature weights

    M::Matrix{Float64} # reference mass matrix
    SS::Vector{Matrix{Float64}} # skew symmetric operator
    MVF::Matrix{Float64} # M^-1 [chiq; chif]^T

    bh::Vector{Vector{Float64}} # boundary integration weights times ref normal
    Ph::Matrix{Float64} # reference projection matrix
    LIFT::Vector{Matrix{Float64}} # lifting operator in each direction

    function RefElemSBP(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes)
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

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Ph = M \ chiq' * Diagonal(wq)

        LIFT = (b -> (M \ chif') * Diagonal(b)).(bh)
        MVF = M \ [chiq ; chif]'

        # We finally build the skew symmetric SBP operator
        SS = Vector{Matrix{Float64}}(undef, dim)
        for idim in 1:dim
            SSi = zeros((Nqnodes + Nfnodes * Nfaces, Nqnodes + Nfnodes * Nfaces))

            Q = Ph' * chiq' * Diagonal(wq) * dchi[idim] * Ph
            EB = (chif * Ph)' * Diagonal(bh[idim])

            SSi[1:Nqnodes, 1:Nqnodes] = Q - Q'
            SSi[1:Nqnodes, Nqnodes+1:end] = EB
            SSi[Nqnodes+1:end, 1:Nqnodes] = -EB'

            SS[idim] = SSi
        end

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
            qnodes,
            fnodes,

            nref,

            chiq,
            chif,
            wq,

            M,
            SS,
            MVF,
            bh,
            Ph,
            LIFT
        )
    end

    function RefElemSBP(bnodestype::AbstractNodes, qnodestype::AbstractNodes, qmnodestype::AbstractNodes, fnodestype::Fnodes)
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

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chimq, wmq, _, _, _ = extract_volume_operators(bnodestype, qmnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Mm = chimq' * Diagonal(wmq) * chimq
        Ph = M \ chiq' * Diagonal(wq)

        LIFT = (b -> (Mm \ chif') * Diagonal(b)).(bh)
        MVF = Mm \ [chiq ; chif]'

        # We finally build the skew symmetric SBP operator
        SS = Vector{Matrix{Float64}}(undef, dim)
        for idim in 1:dim
            SSi = zeros((Nqnodes + Nfnodes * Nfaces, Nqnodes + Nfnodes * Nfaces))

            Q = Ph' * chiq' * Diagonal(wq) * dchi[idim] * Ph
            EB = (chif * Ph)' * Diagonal(bh[idim])

            SSi[1:Nqnodes, 1:Nqnodes] = Q - Q'
            SSi[1:Nqnodes, Nqnodes+1:end] = EB
            SSi[Nqnodes+1:end, 1:Nqnodes] = -EB'

            SS[idim] = SSi
        end

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
            qnodes,
            fnodes,

            nref,

            chiq,
            chif,
            wq,

            M,
            SS,
            MVF,
            bh,
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
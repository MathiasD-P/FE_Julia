# NOTE: CHANGE PACKAGING OF SURFACE TERMS EVENTUALLY.

abstract type AbstractRefElem end

struct RefElemStd <: AbstractRefElem
    shape::String # type of element

    dim::Integer # dimension of element
    Nfaces::Integer # total number of faces
    Nbnodes::Integer # total number of basis nodes
    Nqnodes::Integer # total number of quadrature nodes
    Nqmnodes::Integer # total number of mass quadrature nodes
    Nqsnodes::Integer # total number of source quadrature nodes
    Nfnodes::Integer # number of face nodes on EACH face

    bnodestype::AbstractNodes
    qnodestype::AbstractNodes
    qmnodestype::AbstractNodes
    qsnodestype::AbstractNodes
    fnodestype::Fnodes

    bnodes::Matrix{Float64}
    qnodes::Matrix{Float64}
    qmnodes::Matrix{Float64}
    qsnodes::Matrix{Float64}
    fnodes::Matrix{Float64}

    nref::Vector{Vector{Float64}} # reference normal vectors

    chiq::Matrix{Float64} # volume quadrature Vandermonde matrix
    chiqm::Matrix{Float64} # volume mass quadrature Vandermonde matrix
    chiqs::Matrix{Float64} # volume source quadrature Vandermonde matrix
    chif::Matrix{Float64} # face quadrature Vandermonde matrix
    wq::Vector{Float64} # volume quadrature weights
    wqm::Vector{Float64} # volume mass quadrature weights
    wqs::Vector{Float64} # volume source quadrature weights

    M::Matrix{Float64}# reference mass matrix
    bh::Vector{Vector{Float64}} # boundary integration weights times ref normal
    Dh::Vector{Matrix{Float64}} # reference differentiation matrix in each direction
    Qh::Vector{Matrix{Float64}} # reference stiffness matrix in each direction
    MinvQhT::Vector{Matrix{Float64}} # for weak DG on linear meshes
    Ph::Matrix{Float64} # reference projection matrix
    LIFT::Vector{Matrix{Float64}} # lifting operator in each direction
    Phs::Matrix{Float64} # reference source projection matrix

    function RefElemStd(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes)
        shape = fnodestype.refshape

        dim = bnodestype.dim
        Nfaces = numfaces(fnodestype)
        Nbnodes = numnodes(bnodestype)
        Nqnodes = numnodes(qnodestype)
        Nqmnodes = Nqnodes
        Nqsnodes = Nqnodes
        Nfnodes = numnodes(fnodestype.nodes)

        nref = unitnormals(fnodestype)

        bnodestype = bnodestype
        qnodestype = qnodestype
        qmnodestype = qnodestype
        qsnodestype = qnodestype
        fnodestype = fnodestype

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        qmnodes = qnodes
        chiqm = chiq
        wqm = wq

        qsnodes = qnodes
        chiqs = chiq
        wqs = wq

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Ph = (M \ chiq' * Diagonal(wq))
        Qh = (d -> chiq' * Diagonal(wq) * d).(dchi)
        Dh = (q -> M \ q).(Qh)
        MinvQhT = (q -> M \ q').(Qh)
        LIFT = (b -> (M \ chif') * Diagonal(b)).(bh)

        Phs = Ph

        new(
            shape,

            dim,
            Nfaces,
            Nbnodes,
            Nqnodes,
            Nqmnodes,
            Nqsnodes,
            Nfnodes,

            bnodestype,
            qnodestype,
            qmnodestype,
            qsnodestype,
            fnodestype,

            bnodes,
            qnodes,
            qmnodes,
            qsnodes,
            fnodes,

            nref,

            chiq,
            chiqm,
            chiqs,
            chif,
            wq,
            wqm,
            wqs,

            M,
            bh,
            Dh,
            Qh,
            MinvQhT,
            Ph,
            LIFT,

            Phs
        )
    end

    function RefElemStd(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes, qmnodestype::AbstractNodes, qsnodestype::AbstractNodes)
        shape = fnodestype.refshape

        dim = bnodestype.dim
        Nfaces = numfaces(fnodestype)
        Nbnodes = numnodes(bnodestype)
        Nqnodes = numnodes(qnodestype)
        Nqmnodes = numnodes(qmnodestype)
        Nqsnodes = numnodes(qsnodestype)
        Nfnodes = numnodes(fnodestype.nodes)

        nref = unitnormals(fnodestype)

        bnodestype = bnodestype
        qnodestype = qnodestype
        qmnodestype = qmnodestype
        qsnodestype = qsnodestype
        fnodestype = fnodestype

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chiqm, wqm, _, _, qmnodes = extract_volume_operators(bnodestype, qmnodestype)
        chiqs, wqs, _, _, qsnodes = extract_volume_operators(bnodestype, qsnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Mm = chiqm' * Diagonal(wqm) * chiqm
        Ph = (M \ chiq' * Diagonal(wq))
        Qh = (d -> chiq' * Diagonal(wq) * d).(dchi)
        Dh = (q -> Mm \ q).(Qh)
        MinvQhT = (q -> Mm \ q').(Qh)
        LIFT = (b -> (Mm \ chif') * Diagonal(b)).(bh)

        Phs = Mm \ chiqs' * Diagonal(wqs)

        new(
            shape,

            dim,
            Nfaces,
            Nbnodes,
            Nqnodes,
            Nqmnodes,
            Nqsnodes,
            Nfnodes,

            bnodestype,
            qnodestype,
            qmnodestype,
            qsnodestype,
            fnodestype,

            bnodes,
            qnodes,
            qmnodes,
            qsnodes,
            fnodes,

            nref,

            chiq,
            chiqm,
            chiqs,
            chif,
            wq,
            wqm,
            wqs,

            M,
            bh,
            Dh,
            Qh,
            MinvQhT,
            Ph,
            LIFT,

            Phs
        )
    end

    function RefElemStd(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes, qmnodestype::AbstractNodes)
        return RefElemStd(bnodestype, qnodestype, fnodestype, qmnodestype, qnodestype)
    end
end

struct RefElemSBP <: AbstractRefElem
    shape::String # type of element

    dim::Integer # dimension of element
    Nfaces::Integer # total number of faces
    Nbnodes::Integer # total number of basis nodes
    Nqnodes::Integer # total number of quadrature nodes
    Nqmnodes::Integer # total number of mass quadrature nodes
    Nqsnodes::Integer # total number of source quadrature nodes
    Nfnodes::Integer # number of face nodes on EACH face

    bnodestype::AbstractNodes
    qnodestype::AbstractNodes
    qmnodestype::AbstractNodes
    qsnodestype::AbstractNodes
    fnodestype::Fnodes

    bnodes::Matrix{Float64}
    qnodes::Matrix{Float64}
    qmnodes::Matrix{Float64}
    qsnodes::Matrix{Float64}
    fnodes::Matrix{Float64}

    nref::Vector{Vector{Float64}} # reference normal vectors

    chiq::Matrix{Float64} # volume quadrature Vandermonde matrix
    chiqm::Matrix{Float64} # volume mass quadrature Vandermonde matrix
    chiqs::Matrix{Float64} # volume source quadrature Vandermonde matrix
    chif::Matrix{Float64} # face quadrature Vandermonde matrix
    wq::Vector{Float64} # volume quadrature weights
    wqm::Vector{Float64} # volume mass quadrature weights
    wqs::Vector{Float64} # volume source quadrature weights

    M::Matrix{Float64} # reference mass matrix
    SS::Vector{Matrix{Float64}} # skew symmetric operator
    MVF::Matrix{Float64} # M^-1 [chiq; chif]^T

    bh::Vector{Vector{Float64}} # boundary integration weights times ref normal
    Ph::Matrix{Float64} # reference projection matrix
    LIFT::Vector{Matrix{Float64}} # lifting operator in each direction

    Phs::Matrix{Float64} # reference source projection matrix

    function RefElemSBP(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes)
        shape = fnodestype.refshape

        dim = bnodestype.dim
        Nfaces = numfaces(fnodestype)
        Nbnodes = numnodes(bnodestype)
        Nqnodes = numnodes(qnodestype)
        Nqmnodes = Nqnodes
        Nqsnodes = Nqnodes
        Nfnodes = numnodes(fnodestype.nodes)

        nref = unitnormals(fnodestype)

        bnodestype = bnodestype
        qnodestype = qnodestype
        qmnodestype = qnodestype
        qsnodestype = qnodestype
        fnodestype = fnodestype

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        qmnodes = qnodes
        chiqm = chiq
        wqm = wq

        qsnodes = qnodes
        chiqs = chiq
        wqs = wq

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Ph = M \ chiq' * Diagonal(wq)

        Phs = Ph

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
            Nqmnodes,
            Nqsnodes,
            Nfnodes,

            bnodestype,
            qnodestype,
            qmnodestype,
            qsnodestype,
            fnodestype,

            bnodes,
            qnodes,
            qmnodes,
            qsnodes,
            fnodes,

            nref,

            chiq,
            chiqm,
            chiqs,
            chif,
            wq,
            wqm,
            wqs,

            M,
            SS,
            MVF,
            bh,
            Ph,
            LIFT,

            Phs,
        )
    end

    function RefElemSBP(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes, qmnodestype::AbstractNodes, qsnodestype::AbstractNodes)
        shape = fnodestype.refshape

        dim = bnodestype.dim
        Nfaces = numfaces(fnodestype)
        Nbnodes = numnodes(bnodestype)
        Nqnodes = numnodes(qnodestype)
        Nqmnodes = numnodes(qmnodestype)
        Nqsnodes = numnodes(qsnodestype)
        Nfnodes = numnodes(fnodestype.nodes)

        nref = unitnormals(fnodestype)

        bnodestype = bnodestype
        qnodestype = qnodestype
        qmnodestype = qmnodestype
        qsnodestype = qsnodestype
        fnodestype = fnodestype

        chiq, wq, dchi, bnodes, qnodes = extract_volume_operators(bnodestype, qnodestype)
        chiqm, wqm, _, _, qmnodes = extract_volume_operators(bnodestype, qmnodestype)
        chiqs, wqs, _, _, qsnodes = extract_volume_operators(bnodestype, qsnodestype)
        chif, wf, fnodes = extract_face_operators(bnodestype, fnodestype)

        fnodes = vcat([fnodes[iface] for iface in 1:Nfaces]...)
        chif = vcat([chif[iface] for iface in 1:Nfaces]...) # VERY INEFFICIENT
        bh = [vcat([wf[iface] .* nref[iface][idim] for iface in 1:Nfaces]...) for idim=1:dim]

        M = chiq' * Diagonal(wq) * chiq
        Mm = chiqm' * Diagonal(wqm) * chiqm
        Ph = M \ chiq' * Diagonal(wq)
        Phs = Mm \ chiqs' * Diagonal(wqs)

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
            Nqmnodes,
            Nqsnodes,
            Nfnodes,

            bnodestype,
            qnodestype,
            qmnodestype,
            qsnodestype,
            fnodestype,

            bnodes,
            qnodes,
            qmnodes,
            qsnodes,
            fnodes,

            nref,

            chiq,
            chiqm,
            chiqs,
            chif,
            wq,
            wqm,
            wqs,

            M,
            SS,
            MVF,
            bh,
            Ph,
            LIFT,

            Phs,
        )
    end

    function RefElemSBP(bnodestype::AbstractNodes, qnodestype::AbstractNodes, fnodestype::Fnodes, qmnodestype::AbstractNodes)
        return RefElemSBP(bnodestype, qnodestype, fnodestype, qmnodestype, qnodestype)
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
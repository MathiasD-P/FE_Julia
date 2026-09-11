using LegendrePolynomials
using LinearAlgebra
using FE_Julia

# Expression for Delta_{ijk}(\alpha, \beta) based on paper by Mathias and Siva
# This returns the tensor entry for input (i,j) and k = \tau + 2 - (i + j)
# The input nodes correspond to the choice of quadrature nodes.
function Delta_analytical(nodes::FE_Julia.AbstractNodes, tpflux::TPFlux, i::Integer, j::Integer)
    if nodes.dim != 1
        error("This is a one-dimensional test case!")
    end

    p = nodes.K - 1
    alpha = tpflux.alpha
    beta = tpflux.beta
    
    if nodes.name == "GLL"
        tau = 2 * p -1
        gamma = 2 / (2*p-1) * binomial(4*p,2*p) / binomial(2*p+2,p+1) / binomial(2*p-2,p-1)
    elseif nodes.name == "GL"
        tau = 2 * p + 1
        gamma = -2 / (2*p+3) * binomial(4*p+4, 2*p+2) / binomial(2*p+2,p+1)^2
    end

    delta = binomial(2*(tau+2-(i+j)),tau+2-(i+j)) * binomial(2*i,i) * binomial(2*j,j) / binomial(2*tau+2,tau+1)
    delta *= 0.5 * gamma * (beta * (i+j) - alpha * (tau + 2 - (i+j)))

    return delta
end

# Delta_{i,j,k} is computed directly from our solver
function Delta_empirical(nodes::FE_Julia.AbstractNodes, tpflux::TPFlux, i::Integer, j::Integer, k::Integer)
    if nodes.dim != 1
        error("This is a one-dimensional test case!")
    end

    # Linear mesh with one element
    mesh = FE_Julia.make_interval([-1.0, 1.0], [-1, -1])

    # Build Delta integrals for the split formulation
    dgSBP = DGFluxDiff(1, RefElemSBP(nodes, nodes, make_nodes("interval", "(1)-GLL")), mesh)
    # Fsplit = [FE_Julia.compute_two_pt_flux(Pl(dgSBP.qpts[rowi], i), Pl(dgSBP.qpts[coli], j), tpflux, FE_Julia.Burgers(1)) for rowi in 1:nodes.K, coli in 1:nodes.K]
    samppts = vec([dgSBP.qpts; -1.0; 1.0])
    Fsplit = [0.25 * tpflux.alpha .* (Pl(x,i)*Pl(x,j) + Pl(y,i)*Pl(y,j)) .+ 0.25 * tpflux.beta * (Pl(x,i)*Pl(y,j) .+ Pl(y,i)*Pl(x,j)) for x in samppts, y in samppts]
    ressplit = dgSBP.refelem.M * dgSBP.refelem.MVF * sum(dgSBP.refelem.SS[1] .* Fsplit, dims = 2) .* dgSBP.mesh.detJ[1]

    # Build Delta integrals for the exactly integrated case
    nodesExact = FE_Julia.Tensorprod_nodes("GL", 3 * nodes.K, 1)
    dgExact = DGStd(1, RefElemStd(nodes, nodesExact, make_nodes("interval", "(1)-GLL")), mesh)
    f = 0.5 .* Pl.(vec(dgExact.qpts), i) .* Pl.(vec(dgExact.qpts), j)
    resExact = -dgExact.refelem.Qh[1]' * dgExact.refelem.Ph * f .* dgSBP.mesh.detJ[1]

    # Construct Legendre Vandermonde
    V = [Pl(x, p) for x in vec(dgSBP.bpts), p in 0:(nodes.K-1)]

    # Now, we compute delta
    delta = V' * (ressplit .- resExact)

    return delta[k+1]
end

function test_Delta(nodes::FE_Julia.AbstractNodes, tpflux::TPFlux, i::Integer, j::Integer)
    p = nodes.K - 1

    if nodes.name == "GLL"
        tau = 2 * p -1
    elseif nodes.name == "GL"
        tau = 2 * p + 1
    end

    k = tau + 2 - (i+j)

    if i > p || j > p || k > p
        error("indices are out of range!")
    end
    
    return abs(Delta_analytical(nodes, tpflux, i, j) - Delta_empirical(nodes, tpflux, i, j, k))
end
#####################################################################
# Everything related to the polynomial bases used in HO_FE
# -> Not really optimized, but shouldn't be a too big bottleneck since most arrays are evaluated in
#    the init phase.
#####################################################################

function Lagrange(r::Array, nodes::Array, j::Integer)
    y = ones(Float64, size(r));
    prod = 1.0;
    @inbounds for i in eachindex(nodes)
        if j != i
            y .= (r .- nodes[i]) .* y
            prod *= nodes[j] - nodes[i]
        end
    end
    return y ./ prod
end

function Lagrange_Vandermonde1D(r::Array, nodes::Array)
    r = vec(r)
    nodes = vec(nodes)
    V = Matrix{Float64}(undef, length(r), length(nodes))
    @inbounds @simd for j in eachindex(nodes)
        V[:,j] = Lagrange(r, nodes, j)
    end
    return V
end

function DLagrange(r::Array, nodes::Array, j::Integer)
    dy = 0;
    prod = 1;
    for i in eachindex(nodes)
        if i != j
            prod *= nodes[j] - nodes[i]
            term = ones(size(r))
            for k in eachindex(nodes)
                if (k != i) & (k != j)
                    term = (r .- nodes[k]) .* term
                end
            end
            dy = term .+ dy
        end
    end
    return dy / prod
end

function DLagrange_Vandermonde1D(r::Array, nodes::Array)
    r = vec(r)
    V = Matrix{Float64}(undef, length(r), length(nodes))
    for j in eachindex(nodes)
        V[:,j] = DLagrange(r, nodes, j)
    end
    return V
end

"""
    Vandermonde_tensorprod(V1::Array, V2::array)

Evaluation points are taken to be a "tensor product" of the 1D evaluation points and basis functions are taken to be a flattened tensor product of the bases.
"""
function Vandermonde_tensorprod(V1::Array, V2::Array)
    # M1, N1 = size(V1)
    # M2, N2 = size(V2)

    # V = Matrix{Float64}(undef, M1 * M2, N1 * N2)
    # i = 1
    # @inbounds @simd for i1 in 1:M1
    #     @inbounds @simd for i2 in 1:M2
    #         j = 1
    #         @inbounds @simd for j1 in 1:N1
    #             @inbounds @simd for j2 in 1:N2
    #                 V[i,j] = V1[i1,j1] * V2[i2, j2]
    #                 j += 1
    #             end
    #         end
    #         i += 1
    #     end
    # end
    # return V
    return kron(V1, V2)
end

"""
    point_tensorprod(p1::Matrix, p2::Matrix)

Tensor product of points coordinate which is consistent with kron. Returns a matrix, rows contain the points, columns contain the value of
points along each dimension.
"""
function point_tensorprod(p1::Matrix, p2::Matrix)
    M1 = size(p1,1)
    M2 = size(p2,1)

    return [repeat(p1, inner=(M2,1)) repeat(p2, outer=(M1,1))]

    # p = Matrix{Float64}(undef, M1 * M2, N1 + N2)
    # i = 1
    # @inbounds for i1 in 1:M1
    #     @inbounds for i2 in 1:M2
    #         p[i,:] = [p1[i1, :] p2[i2, :]]
    #         i += 1
    #     end
    # end
    # return p
end

function row_kron_matmul(A,B,x)
    M = size(A,2)
    N = size(B,2)
    row = size(A,1)

    x = reshape(x,N,M)
    y = Vector{Float64}(undef, row)

    @inbounds for i in 1:row
        s = 0.0
        for j in 1:M
            s1 = 0.0
            @simd for k in 1:N
                s1 += B[i,k] * x[k,j]
            end
            s1 *= A[i,j]
            s += s1
        end
        y[i] = s
    end
    return y
end
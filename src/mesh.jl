#####################################################################
# Create, manipulate and store meshes for High-order FE purposes
#####################################################################

abstract type AbstractMesh end

"""
    LMesh: Main object used to store and manipulate infromation about the mesh (for 1 and 2D)
    - elements: (Nel, Nvertex) Matrix of vectors storing the the vertices of each element.
                Induces a global numbering for the elements and a local ordering for each vertex with respect
                to the reference element. Vertices must be numbered counter-clockwise.
    - connectivity: (Nel, Nfaces) Matrix storing the "connectivity graph" of the mesh. Entries are integer vectors of
                    the form [ielem, reffaceindex].
"""
mutable struct LMesh <: AbstractMesh
    dim::Integer # dimension of the space discretized
    Nvert::Integer # Number of vertices per element
    Nface::Integer # Number of faces per element
    Nel::Integer # Number of elements
    Ntags::Integer # Number of BC tags

    elements::AbstractMatrix
    connectivity::AbstractMatrix
    BCtags::Vector{Integer} # MUST BE NEGATIVE
    
    J::Vector{Matrix} # Jacobian for the forward transformation xi -> x
    detJ::Vector # determinant of the Jacobian(s)
    CT::Vector{Matrix} # cofactor matrix transposed
end

"""
    mapping(mesh::LMesh, refpoints::Matrix, ielem::Integer)
    maps the array of ref points (Npoints, dim) to a new array of the same size.
"""
function mapping(mesh::LMesh, refpoints::Matrix, ielem::Integer)
    N = size(refpoints,1)
    return refpoints * mesh.J[ielem]' + repeat(mesh.elements[ielem, 1], 1, N)'
end

function invmapping(mesh::LMesh, points::Matrix, ielem::Integer) # Should store inverse somewhere instead
    N = size(points,1)
    return (points - repeat(mesh.elements[ielem, 1], 1, N)') * (mesh.CT[ielem])' ./ mesh.detJ[ielem]
end

function make_rectangle_quad(Nx::Integer, Ny::Integer, BCtags = [0, 0, 0, 0], Lx=1, Ly=1, origin::Vector = [0,0])
    hx = Lx / Nx
    hy = Ly / Ny
    Nel = Nx * Ny

    X = collect(range(origin[1], origin[1] + Lx - hx, Nx))
    Y = collect(range(origin[2], origin[2] + Ly - hy, Ny))
    
    elements = Matrix{Vector{Float64}}(undef, Nel, 4)
    connectivity = Matrix{Vector{Int64}}(undef, Nel, 4)

    ielem = 0
    for i in eachindex(X)
        for j in eachindex(Y)
            ielem += 1;
            elements[ielem,1] = [X[i],Y[j]]
            elements[ielem,2] = [X[i]+hx,Y[j]]
            elements[ielem,3] = [X[i]+hx,Y[j]+hy]
            elements[ielem,4] = [X[i],Y[j]+hy]
            
            if j == 1
                if BCtags[1] == 0
                    connectivity[ielem, 1] = [ielem+Ny-1, 3]
                else
                    connectivity[ielem, 1] = [0,BCtags[1]]
                end
            else
                connectivity[ielem, 1] = [ielem-1, 3]
            end
            if j == Ny
                if BCtags[2] == 0
                    connectivity[ielem, 3] = [ielem-Ny+1, 1]
                else
                    connectivity[ielem, 3] = [0, BCtags[2]]
                end
            else
                connectivity[ielem, 3] = [ielem+1, 1]
            end

            if i == 1
                if BCtags[3] == 0
                    connectivity[ielem, 4] = [j + (Nx-1)*Ny, 2]
                else
                    connectivity[ielem, 4] = [0,BCtags[3]]
                end
            else
                connectivity[ielem, 4] = [ielem - Ny, 2]
            end
            if i == Nx
                if BCtags[4] == 0
                    connectivity[ielem, 2] = [j, 4]
                else
                    connectivity[ielem, 2] = [0,BCtags[4]]
                end
            else
                connectivity[ielem, 2] = [ielem + Ny, 4]
            end
        end
    end

    J = [hx 0 ; 0 hy]
    J = [J for _ in 1:Nel]
    detJ = [hx * hy for _ in 1:Nel]
    CT = [hy 0 ; 0 hx]
    CT = [CT for _ in 1:Nel]

    dim = 2
    Nvert = 4
    Nface = 4

    BCtags = unique(BCtags)
    Ntags = sum(BCtags .!= 0)

    return LMesh(dim, Nvert, Nface, Nel, Ntags, elements, connectivity, BCtags, J, detJ, CT)
end

function make_interval(X::Vector, BCtags=[0, 0])
    Nel = length(X)-1
    h = X[2:end] - X[1:end-1]

    elements = Matrix{Vector{Float64}}(undef, Nel, 2)
    connectivity = Matrix{Vector{Int64}}(undef, Nel, 2)

    for ielem in 1:Nel
        elements[ielem,1] = [X[ielem]]
        elements[ielem,2] = [X[ielem+1]]
        
        if ielem == 1
            if BCtags[1] == 0
                connectivity[ielem,1] = [Nel, 2]
            else
                connectivity[ielem,1] = [0, BCtags[1]]
            end
        else
            connectivity[ielem,1] = [ielem-1, 2]
        end
        if ielem == Nel
            if BCtags[2] == 0
                connectivity[ielem,2] = [1, 1]
            else
                connectivity[ielem,2] = [0, BCtags[2]]
            end
        else
            connectivity[ielem,2] = [ielem+1, 1]
        end
    end

    J = [reshape([h[i]], (1,1)) for i in 1:Nel]
    detJ = vec(h)

    CT = reshape([1.0], (1,1))
    CT = [CT for _ in 1:Nel]

    dim = 1
    Nvert = 2
    Nface = 2
    
    BCtags = unique(BCtags)
    Ntags = sum(BCtags .!= 0)

    return LMesh(dim, Nvert, Nface, Nel, Ntags, elements, connectivity, BCtags, J, detJ, CT)
end
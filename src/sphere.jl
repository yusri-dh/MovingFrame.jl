###
# using LinearAlgebra:sqrt
# using GeometryTypes
###
###
"""
Return vertex coordinates fixed to the unit sphere
"""
function vertex(x::T, y::T, z::T,scale) where T
    length = sqrt(x^2 + y^2 + z^2)
    return [(i * scale) / length for i in (x,y,z)]
end

function vertex(a::Array{T,1},scale) where T
    @assert length(a) == 3
    x,y,z = a[1],a[2],a[3]
    return vertex(x, y, z,scale)
end

"""
Find a middle point and project to the unit sphere
"""
function middle_point!(point_1,point_2,middle_point_cache,verts,scale)
    # We check if we have already cut this edge first
    # to avoid duplicated verts
    smaller_index = min(point_1, point_2)
    greater_index = max(point_1, point_2)
    key = "$smaller_index, $greater_index"
    if haskey(middle_point_cache,key)
        return middle_point_cache[key]
    end
    vert_1 = verts[point_1]
    vert_2 = verts[point_2]
    middle = [sum(i)/2 for i in zip(vert_1, vert_2)]
    append!(verts,[vertex(middle,scale)])
    index = length(verts)
    middle_point_cache[key] = index
    return index
end # function

function middle_point!(point_1::T,point_2::T,middle_point_cache,verts) where {T<:OffsetInteger}
    # to avoid duplicated verts
    smaller_index = min(point_1.i, point_2.i)
    greater_index = max(point_1.i, point_2.i)
    key = "$smaller_index, $greater_index"
    if haskey(middle_point_cache,key)
        return middle_point_cache[key]
    end
    vert_1 = verts[point_1]
    vert_2 = verts[point_2]
    middle = [sum(i)/2 for i in zip(vert_1, vert_2)]
    append!(verts,[vertex(middle,scale)])
    index = length(verts)
    middle_point_cache[key] = index
    return index
end # function


"""
Create icosphere
"""
function icosphere(subdiv=3,scale = 1.0)
    PHI = (1. + sqrt(5.)) / 2.
    verts = [(vertex(-1., PHI, 0.,scale)),
        (vertex(1., PHI, 0.,scale)),
        (vertex(-1., -PHI, 0.,scale)),
        (vertex(1., -PHI, 0.,scale)),
        (vertex(0., -1., PHI,scale)),
        (vertex(0., 1., PHI,scale)),
        (vertex(0., -1., -PHI,scale)),
        (vertex(0., 1., -PHI,scale)),
        (vertex(PHI, 0., -1.,scale)),
        (vertex(PHI, 0., 1.,scale)),
        (vertex(-PHI, 0., -1.,scale)),
        (vertex(-PHI, 0., 1.,scale))]
    faces =[[1, 12, 6], # 5 faces around point 0
            [1, 6, 2],
            [1, 2, 8],
            [1, 8, 11],
            [1, 11, 12],
            [2, 6, 10], # Adjacent faces
            [6, 12, 5],
            [12, 11, 3],
            [11, 8, 7],
            [8, 2, 9],
            [4, 10, 5],  # 5 faces around 3
            [4, 5, 3],
            [4, 3, 7],
            [4, 7, 9],
            [4, 9, 10],
            [5, 10, 6],# Adjacent faces
            [3, 5, 12],
            [7, 3, 11],
            [9, 7, 8],
            [10, 9, 2]]
    
    middle_point_cache = Dict{String,Int}()
    for i in 1:subdiv
        faces_subdiv = Array{Int64,1}[]
        for tri in faces
            v1 = middle_point!(tri[1], tri[2],middle_point_cache,verts,scale)
            v2 = middle_point!(tri[2], tri[3],middle_point_cache,verts,scale)
            v3 = middle_point!(tri[3], tri[1],middle_point_cache,verts,scale)
            append!(faces_subdiv,[[tri[1], v1, v3]])
            append!(faces_subdiv,[[tri[2], v2, v1]])
            append!(faces_subdiv,[[tri[3], v3, v2]])
            append!(faces_subdiv,[[v1, v2, v3]])
            faces = faces_subdiv
        end
    end
    v = array_to_vertex(verts)
    f = array_to_faces(faces)
    sphere = GLNormalMesh(v,f)
    return sphere
end

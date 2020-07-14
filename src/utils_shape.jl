# using LinearAlgebra
function cart2sph(x::T, y::T, z::T) where {T<:Real}
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    pol = atan(hxy, z)
    az = atan(y, x)
    return [Float64(pol), Float64(az), Float64(r)]
end

function cart2sph(xyz::Point)
    x, y, z = xyz[1], xyz[2], xyz[3]
    return cart2sph(x, y, z)
end

function sph2cart(pol::T, az::T, r::T) where {T<:Real}
    x = r * sin(pol) * cos(az)
    y = r * sin(pol) * sin(az)
    z = r * cos(pol)
    return x, y, z
end
"""
Change array of vertex coordinates into array of vertex type
"""
function array_to_vertex(a::AbstractArray)
    return [Point3f0(vertex) for vertex in a]
end

"""
Change array of face into array of face type
"""
function array_to_faces(a::AbstractArray)
    return [TriangleFace{OffsetInteger{-1,UInt32}}(face) for face in a]
end

"""
Change array of array of into matrix
"""
function array_of_array_to_matrix(a::AbstractArray)
    return [a[i][j] for i = 1:length(a), j = 1:length(a[1])]
end
"""
Change matrix to array of array
"""
function matrix_to_array_of_array(mat::AbstractMatrix)
    return [mat[i, :] for i = 1:size(mat, 1)]
end # function]

"""
Change matrix to vertices
"""
function matrix_to_vertex(mat::AbstractMatrix)
    return [Point3f0(mat[i, :]) for i = 1:size(mat, 1)]
end # function


"""
Centering the coordinates of mesh vertices to their mean
"""
function centroid(v::AbstractMatrix)
    return sum(v, dims = 1) / size(v, 1)
end
function centroid(verts::Vector{T1}) where {T1}
    return reduce(+, verts) / length(verts)
end

function mesh_centering(verts::Vector{T1}, faces::Vector{T2}) where {T1,T2}
    mean_coord = reduce(+, verts) / length(verts)
    new_verts = similar(verts)
    for i in eachindex(verts)
        new_verts[i] = verts[i] - mean_coord
    end
    return new_verts, faces
end
function mesh_centering(mesh::Mesh)
    verts = decompose(Point3f0, mesh)
    triangles = decompose(TriangleFace{Int}, mesh)
    verts, triangles = mesh_centering(verts, triangles)
    return Mesh(verts, faces(mesh))
end



function signed_volume_of_triangle(
    p1::T,
    p2::T,
    p3::T,
) where {T<:AbstractVector}
    return dot(p1, cross(p2, p3)) / 6.0
end


"""
These functions return the volume of mesh object
"""
function mesh_volume(
    verts::Array{Point{3,T1}},
    faces::Array{TriangleFace{T2}},
) where {T1<:Real,T2<:Int}
    total_vol = 0.0
    for face in faces
        vol = signed_volume_of_triangle(
            verts[face[1]],
            verts[face[2]],
            verts[face[3]],
        )
        total_vol += vol
    end
    return abs(total_vol)
end

function mesh_volume(mesh::Mesh)
    verts = decompose(Point3f0, mesh)
    faces = decompose(TriangleFace{Int}, mesh)

    return mesh_volume(verts, faces)
end
"""
These functions return the copy of mesh object with volume one
"""
function volume_normalizing(
    verts::Array{Point{3,T1}},
    faces::Array{TriangleFace{T2}},
) where {T1<:Real,T2<:Int}
    vol = mesh_volume(verts, faces)
    cbrt_vol = cbrt(vol)
    new_verts = similar(verts)
    for i in eachindex(verts)
        new_verts[i] = verts[i] / cbrt_vol
    end
    return new_verts, faces
end

function volume_normalizing2(
    verts::Array{Point{3,T1}},
    faces::Array{TriangleFace{T2}},
) where {T1<:Real,T2<:Int}
    vol = mesh_volume(verts, faces)
    cbrt_vol = cbrt(vol)
    new_verts = similar(verts)
    for i in eachindex(verts)
        new_verts[i] = verts[i] / cbrt_vol
    end
    return new_verts, faces
end

function volume_normalizing(mesh::Mesh)
    verts = decompose(Point3f0, mesh)
    triangles = decompose(TriangleFace{Int}, mesh)
    verts, triangles = volume_normalizing(verts, triangles)
    return Mesh(verts, faces(mesh))
end
"""
This function return the mesh object from verts and faces
"""
function create_mesh(verts, faces)
    faces = convert(Vector{TriangleFace{OffsetInteger{-1,UInt32}}}, faces)
    return Mesh(verts, faces)
end
"""
This function return the angle at vertex v1 from triangle(v1,v2,v3)
"""
function calculate_angle(v1, v2, v3)
    v12 = v2 - v1
    v13 = v3 - v1
    v12 = v12 / norm(v12)
    v13 = v13 / norm(v13)
    cos_val = dot(v12, v13)
    angle = acos(cos_val)
    return angle
end
# """
# This function return the internal angles of vertices
# """
# function internal_angles(verts,faces)
#     nrow,ncol = length(faces[1]),length(faces)
#     angles = Matrix{Float64}(undef,nrow,ncol)
#     for j in 1:ncol
#         f = faces[j]
#         for  i in 1:nrow
#             v1 = verts[f[i]]
#             idx2 = mod(i+1,nrow) == 0 ? 3 : mod(i+1,nrow)
#             idx3 = mod(i-1,nrow) == 0 ? 3 : mod(i-1,nrow)
#             # @show i
#             # @show idx3
#             # @show f
#             v2 = verts[f[idx2]]
#             v3 = verts[f[idx3]]
#
#             angles[i,j] = calculate_angle(v1,v2,v3)
#         end
#     end
#     return angles
# end
"""
This function return the internal angles of vertices
"""
function internal_angles_intrinsic(l)
    s23 = l[:, 1]
    s31 = l[:, 2]
    s12 = l[:, 3]

    a23 = acos.((s12 .^ 2 + s31 .^ 2 - s23 .^ 2) ./ (2 * s12 .* s31))
    a31 = acos.((s23 .^ 2 + s12 .^ 2 - s31 .^ 2) ./ (2 * s23 .* s12))
    a12 = acos.((s31 .^ 2 + s23 .^ 2 - s12 .^ 2) ./ (2 * s31 .* s23))

    angles = [a23 a31 a12]
    return angles
end
function internal_angles(verts, faces)
    f1 = [face[1] for face in faces]
    f2 = [face[2] for face in faces]
    f3 = [face[3] for face in faces]

    s12 = norm.(verts[f2] - verts[f1])
    s13 = norm.(verts[f3] - verts[f1])
    s23 = norm.(verts[f3] - verts[f2])

    l = [s23 s13 s12]

    A = internal_angles_intrinsic(l)
end

"""
This function return the area of triangle made from vectors v1,v2, and v3)
"""
function calculate_area(v1, v2, v3)
    v12 = v2 - v1
    v13 = v3 - v1
    area = norm(cross(v12, v13)) / 2.0
    return area
end
"""
This function return the surface area)
"""
function surface_area(v::AbstractVector, f::AbstractVector)
    nv = length(v)
    f1 = [face[1] for face in f]
    f2 = [face[2] for face in f]
    f3 = [face[3] for face in f]

    l1 = norm.(v[f2] - v[f3])
    l2 = norm.(v[f3] - v[f1])
    l3 = norm.(v[f1] - v[f2])
    s = (l1 + l2 + l3) / 2.0
    area = sqrt.(s .* (s - l1) .* (s - l2) .* (s - l3))
    return sum(area)
end
function surface_area(v::Matrix{<:AbstractFloat}, f::Matrix{<:Signed})
    f1 = @view f[:, 1]
    f2 = @view f[:, 2]
    f3 = @view f[:, 3]

    l1 = sqrt.(sum((v[f2, :] - v[f3, :]) .^ 2, dims = 2))
    l2 = sqrt.(sum((v[f3, :] - v[f1, :]) .^ 2, dims = 2))
    l3 = sqrt.(sum((v[f1, :] - v[f2, :]) .^ 2, dims = 2))
    s = (l1 + l2 + l3) / 2.0
    # @show s
    area = sqrt.(s .* (s - l1) .* (s - l2) .* (s - l3))
    return sum(area)
end
function surface_area2(v::AbstractVector, f::AbstractVector)
    f1 = [face[1] for face in f]
    f2 = [face[2] for face in f]
    f3 = [face[3] for face in f]
    #F x 3 matrices of triangle edge vectors, named after opposite vertices
    v1 = [v[f3[i]] - v[f2[i]] for i in eachindex(f)]
    v2 = [v[f1[i]] - v[f3[i]] for i in eachindex(f)]
    # v3 = [v[f2[i]] - v[f1[i]] for i in eachindex(f)];
    n = cross.(v1, v2)
    area = norm.(n) / 2.0
    return sum(area)
end
function surface_area2(v::Matrix{<:AbstractFloat}, f::Matrix{<:Signed})
    f1 = @view f[:, 1]
    f2 = @view f[:, 2]
    f3 = @view f[:, 3]
    #F x 3 matrices of triangle edge vectors, named after opposite vertices
    v1 = v[f3, :] - v[f2, :]
    v2 = v[f1, :] - v[f3, :]
    # v3 = v[f2,:] - v[f1,:];
    n = [cross(v1[i, :], v2[i, :]) for i = 1:size(v1, 1)]
    area = norm.(n) / 2.0
    return sum(area)
end
"""
This function return the surface area)
"""
function double_surface_area(v, f)
    return surface_area(v, f) * 2.0
end

function mesh_reorientation(mesh::Mesh, T, N, B)
    verts = decompose(Point3f0, mesh)
    mat_verts = array_of_array_to_matrix(verts)
    new_mat_verts = [T N B] \ mat_verts'
    new_mat_verts = new_mat_verts'
    new_verts = matrix_to_vertex(new_mat_verts)
    return Mesh(new_verts, faces(mesh))
end
function mesh_reorientation(mat_verts::AbstractMatrix, T, N, B)
    new_mat_verts = [T N B] \ mat_verts'
    return new_mat_verts'
end

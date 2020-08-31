function mass_matrix_barycentric(v, f)
    f1 = [face[1] for face in f]
    f2 = [face[2] for face in f]
    f3 = [face[3] for face in f]
    # F x 3 matrices of triangle edge vectors, named after opposite vertices
    v1 = [v[f3[i]] - v[f2[i]] for i in eachindex(f)];
    v2 = [v[f1[i]] - v[f3[i]] for i in eachindex(f)];
    v3 = [v[f2[i]] - v[f1[i]] for i in eachindex(f)];
    n = cross.(v1, v2)
    area = norm.(n) / 2
    i = [f1;f2;f3]
    j = [f1;f2;f3]


    diag_v = area / 3
    val = [diag_v;diag_v;diag_v];
    M = sparse(i, j, val, size(v, 1), size(v, 1));
    return M
end

function mass_matrix_barycentric(v::Matrix{<: AbstractFloat}, f::Matrix{<: Signed})
    f1 = @view f[:,1]
    f2 = @view f[:,2]
    f3 = @view f[:,3]

    l1 = sqrt.(sum((v[f2,:] - v[f3,:]).^2, dims=2))
    l2 = sqrt.(sum((v[f3,:] - v[f1,:]).^2, dims=2))
    l3 = sqrt.(sum((v[f1,:] - v[f2,:]).^2, dims=2))
    s = (l1 + l2 + l3) / 2.0
    area = sqrt.(s .* (s - l1) .* (s - l2) .* (s - l3))
    i = [f1;f2;f3]
    j = [f1;f2;f3]


    diag_v = reshape(area / 3, :)
    val = [diag_v;diag_v;diag_v];
    M = sparse(i, j, val, size(v, 1), size(v, 1));
    return M
end

function mass_matrix_barycentric2(v::Matrix{<: AbstractFloat}, f::Matrix{<: Signed})
    f1 = @view f[:,1]
    f2 = @view f[:,2]
    f3 = @view f[:,3]

    v1 = v[f3,:] - v[f2,:];
    v2 = v[f1,:] - v[f3,:];
    v3 = v[f2,:] - v[f1,:];
    n = [cross(v1[i,:], v2[i,:]) for i in 1:size(v1, 1)]
    area = norm.(n) / 2
    i = [f1;f2;f3]
    j = [f1;f2;f3]


    diag_v = area / 3
    val = [diag_v;diag_v;diag_v];
    M = sparse(i, j, val, size(v, 1), size(v, 1));
    return M
end
"""
This function return the cotangent weighted laplacian matrix.
"""
function cotangent_laplacian(v, f)
    nv = length(v)
    f1 = [face[1] for face in f]
    f2 = [face[2] for face in f]
    f3 = [face[3] for face in f]

    l1 = norm.(v[f2] - v[f3])
    l2 = norm.(v[f3] - v[f1])
    l3 = norm.(v[f1] - v[f2])
    s = (l1 + l2 + l3) / 2.0
    area = sqrt.(s .* (s - l1) .* (s - l2) .* (s - l3))
    cot12 = (l1.^2 + l2.^2 - l3.^2) ./ area / 2
    cot23 = (l2.^2 + l3.^2 - l1.^2) ./ area / 2
    cot31 = (l1.^2 + l3.^2 - l2.^2) ./ area / 2
    diag1 = -cot12 - cot31
    diag2 = -cot12 - cot23
    diag3 = -cot31 - cot23

    II = [f1; f2; f2; f3; f3; f1; f1; f2; f3]
    JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3]
    V = [cot12; cot12; cot23; cot23; cot31; cot31; diag1; diag2; diag3]
    L = sparse(II, JJ, V, nv, nv)
    return L
end
function cotangent_laplacian(v::Matrix{<:AbstractFloat}, f::Matrix{<:Signed})
    nv = size(v, 1)
    f1 = @view f[:,1]
    f2 = @view f[:,2]
    f3 = @view f[:,3]

    l1 = reshape(sqrt.(sum((v[f2,:] - v[f3,:]).^2, dims=2)), :)
    l2 = reshape(sqrt.(sum((v[f3,:] - v[f1,:]).^2, dims=2)), :)
    l3 = reshape(sqrt.(sum((v[f1,:] - v[f2,:]).^2, dims=2)), :)
    s = (l1 + l2 + l3) / 2.0
    area = sqrt.(s .* (s - l1) .* (s - l2) .* (s - l3))
    cot12 = (l1.^2 + l2.^2 - l3.^2) ./ area / 2
    cot23 = (l2.^2 + l3.^2 - l1.^2) ./ area / 2
    cot31 = (l1.^2 + l3.^2 - l2.^2) ./ area / 2
    diag1 = -cot12 - cot31
    diag2 = -cot12 - cot23
    diag3 = -cot31 - cot23

    II = [f1; f2; f2; f3; f3; f1; f1; f2; f3]
    JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3]
    V = [cot12; cot12; cot23; cot23; cot31; cot31; diag1; diag2; diag3]
    L = sparse(II, JJ, V, nv, nv)
    return L
end

function vertex_adjacent_area(verts, faces)
    areas = zeros(Float64, length(verts))
    for j in 1:length(faces)
        f = faces[j]
        area = calculate_area(verts[f[1]], verts[f[2]], verts[f[3]])
        for i in 1:3
            areas[f[i]] += area / 3.0
        end
    end
    return areas
end
"""
This function return the gaussian curvature of the vertex
"""
function gaussian_curvature(verts, faces)

    int_angles = internal_angles(verts, faces)
    gauss_K = repeat([2 * π], length(verts))
    for j in 1:length(faces)
        for i in 1:3
            vert = faces[j][i]
            gauss_K[vert] -= int_angles[j,i]
        end
    end
    # return gauss_K ./ vertex_adjacent_area(verts,faces)
    return gauss_K
end
function gaussian_curvature(mesh::Mesh)
    verts = decompose(Point3f0, mesh)
    faces = decompose(TriangleFace{Int}, mesh)
    return gaussian_curvature(verts, faces)
end
"""
This function return the mean curvature of the vertex
"""
function mean_curvature(verts, faces)
    l = cotangent_laplacian(verts, faces)
    m = mass_matrix_barycentric(verts, faces)
    v = MovingFrame.array_of_array_to_matrix(verts)
    minv = spdiagm(0 => (1.0 ./ diag(m)))

    hn = -minv * (l * v )
    h = sqrt.(sum(hn.^2, dims=2))
    # return gauss_K ./ vertex_adjacent_area(verts,faces)
    return reshape(h, :)
end
function mean_curvature(mesh::Mesh)
    verts = decompose(Point3f0, mesh)
    faces = decompose(TriangleFace{Int}, mesh)
    return mean_curvature(verts, faces)
end
function mean_curvature_flow(verts, faces; n_step=10, delta=0.001,min_diff=1e-13)
    new_verts, _ = volume_normalizing(mesh_centering(verts, faces)...)
    v = MovingFrame.array_of_array_to_matrix(new_verts)
    f = MovingFrame.array_of_array_to_matrix(faces)
    new_v = v .- centroid(v)
    new_v = new_v / sqrt(surface_area2(new_v, f))
    l = cotangent_laplacian(new_v, f)

    for i in 1:n_step
        old_v = copy(new_v)
        m = mass_matrix_barycentric2(new_v, f)
        s = (m - delta * l)
        b = m * new_v
        new_v = s \ b
        new_v = new_v .- centroid(new_v)
        new_v = new_v / sqrt(surface_area2(new_v, f))
        diff = tr((new_v - old_v)' * mass_matrix_barycentric2(new_v, f) * (new_v - old_v))
        # @show diff
        if diff <= min_diff
            println("Converged")
            break
        end
    end
    new_verts::typeof(verts) = MovingFrame.matrix_to_vertex(new_v)
    return new_verts
end
function spherization(verts, faces;n_step=10, delta=0.001,min_diff=1e-13)
    new_verts =  mean_curvature_flow(verts, faces; n_step=n_step, delta=delta,min_diff=min_diff)
    return new_verts, faces
end
function spherization(mesh::Mesh;n_step=10, delta=0.001,min_diff=1e-13)
    verts = decompose(Point3f0, mesh)
    triangles = decompose(TriangleFace{Int}, mesh)
    new_verts =  mean_curvature_flow(verts, triangles; n_step=n_step, delta=delta,min_diff=min_diff)
    return Mesh(new_verts, faces(mesh))
end

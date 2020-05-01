function spharm_descriptor(
    verts::Array{Point{3,T1}},
    faces::Array{Face{3,T2}},
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13
)  where {T1<:Real,T2<:Int}
    if normalized
        verts, faces = volume_normalizing(verts, faces)
    end
    if centered
        verts, faces = mesh_centering(verts, faces)
    end
    verts_sphere, _ = spherization(
        verts,
        faces,
        n_step,
        delta,
        min_diff = min_diff,
    )
    verts_sphere_spherical_coords = cart2sph.(verts_sphere)
    for i = 1:length(verts_sphere_spherical_coords)
        if verts_sphere_spherical_coords[i][2] < 0
            verts_sphere_spherical_coords[i][2] += 2pi
        end
    end
    x = [v[1] for v in verts]
    y = [v[2] for v in verts]
    z = [v[3] for v in verts]
    pol_angle = [v_sp[1] for v_sp in verts_sphere_spherical_coords]
    az_angle = [v_sp[2] for v_sp in verts_sphere_spherical_coords]
    c_x = spharm_coefs(l_max::Int, x, pol_angle, az_angle)
    c_y = spharm_coefs(l_max::Int, y, pol_angle, az_angle)
    c_z = spharm_coefs(l_max::Int, z, pol_angle, az_angle)
    return (Cx = c_x, Cy = c_y, Cz = c_z)
end
function spharm_descriptor(
    mesh::Mesh;
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13,
)
    verts = decompose(Point3f0, mesh)
    faces = decompose(Face{3,Int}, mesh)
    return spharm_descriptor(
        verts,
        faces;
        centered = centered,
        normalized = normalized,
        l_max = l_max,
        n_step = n_step,
        delta = delta,
        min_diff = min_diff,
    )
end
function spharm_descriptor(
    verts::Vector{Point{3,T1}},
    verts_sphere::Vector{Point{3,T1}},
    faces::Vector{Face{3,T2}};
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13
) where {T1<:Real,T2<:Int}
    if normalized
        verts, faces = volume_normalizing(verts, faces)
    end
    if centered
        verts, faces = mesh_centering(verts, faces)
    end
    verts_sphere_spherical_coords = cart2sph.(verts_sphere)
    for i = 1:length(verts_sphere_spherical_coords)
        if verts_sphere_spherical_coords[i][2] < 0
            verts_sphere_spherical_coords[i][2] += 2pi
        end
    end
    x = [v[1] for v in verts]
    y = [v[2] for v in verts]
    z = [v[3] for v in verts]
    pol_angle = [v_sp[1] for v_sp in verts_sphere_spherical_coords]
    az_angle = [v_sp[2] for v_sp in verts_sphere_spherical_coords]
    c_x = spharm_coefs(l_max::Int, x, pol_angle, az_angle)
    c_y = spharm_coefs(l_max::Int, y, pol_angle, az_angle)
    c_z = spharm_coefs(l_max::Int, z, pol_angle, az_angle)
    return (Cx = c_x, Cy = c_y, Cz = c_z)
end
function spharm_descriptor(
    mesh::Mesh,
    verts_sphere::Vector{Point{3,T1}};
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13
) where{T1<:Real}
    verts = decompose(Point3f0, mesh)
    faces = decompose(Face{3,Int}, mesh)
    return spharm_descriptor(
        verts,
        verts_sphere,
        faces,
        centered = centered,
        normalized = normalized,
        l_max = l_max,
        n_step = n_step,
        delta = delta,
        min_diff = min_diff,
    )
end
function spharm_descriptor_v(
    verts::Vector{Point{3,T1}},
    faces::Vector{Face{3,T2}},
    val::Vector{<:AbstractFloat};
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13,
) where {T1<:Real,T2<:Int}
    if normalized
        verts, faces = volume_normalizing(verts, faces)
    end
    if centered
        verts, faces = mesh_centering(verts, faces)
    end
    verts_sphere, _ = spherization(
        verts,
        faces,
        n_step,
        delta,
        min_diff = min_diff,
    )
    verts_sphere_spherical_coords = cart2sph.(verts_sphere)
    for i = 1:length(verts_sphere_spherical_coords)
        if verts_sphere_spherical_coords[i][2] < 0
            verts_sphere_spherical_coords[i][2] += 2pi
        end
    end
    pol_angle = [v_sp[1] for v_sp in verts_sphere_spherical_coords]
    az_angle = [v_sp[2] for v_sp in verts_sphere_spherical_coords]
    c_v = spharm_coefs(l_max::Int, val, pol_angle, az_angle)
    return c_v
end
function spharm_descriptor_v(
    mesh::Mesh,
    val::Vector{<:AbstractFloat};
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13,
)
    verts = decompose(Point3f0, mesh)
    faces = decompose(Face{3,Int}, mesh)
    return spharm_descriptor_v(
        verts,
        faces,
        val::Vector{<:AbstractFloat};
        centered = true,
        normalized = true,
        l_max = 6,
        n_step = 10,
        delta = 0.001,
        min_diff = 1e-13,
    )
end

function spharm_descriptor_v(
    verts::Vector{Point{3,T1}},
    verts_sphere::Vector{Point{3,T1}},
    faces::Vector{Face{3,T2}},
    val::Vector{<:AbstractFloat},;
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13,
) where {T1<:Real,T2<:Int}
    if normalized
        verts, faces = volume_normalizing(verts, faces)
    end
    if centered
        verts, faces = mesh_centering(verts, faces)
    end
    verts_sphere_spherical_coords = cart2sph.(verts_sphere)
    for i = 1:length(verts_sphere_spherical_coords)
        if verts_sphere_spherical_coords[i][2] < 0
            verts_sphere_spherical_coords[i][2] += 2pi
        end
    end
    pol_angle = [v_sp[1] for v_sp in verts_sphere_spherical_coords]
    az_angle = [v_sp[2] for v_sp in verts_sphere_spherical_coords]
    c_v = spharm_coefs(l_max::Int, val, pol_angle, az_angle)
    return c_v
end
function spharm_descriptor_v(
    mesh::Mesh,
    verts_sphere::Vector{Point{3,T1}},
    val::Vector{<:AbstractFloat};
    centered = true,
    normalized = true,
    l_max = 6,
    n_step = 10,
    delta = 0.001,
    min_diff = 1e-13,
) where {T1<:Real}
    verts = decompose(Point3f0, mesh)
    faces = decompose(Face{3,Int}, mesh)
    return spharm_descriptor_v(
        verts,
        verts_sphere,
        faces,
        val::Vector{<:AbstractFloat};
        centered = true,
        normalized = true,
        l_max = 6,
        n_step = 10,
        delta = 0.001,
        min_diff = 1e-13,
    )
end



function reconstruct_from_spharm_coefs(
    coefs::AbstractVector,
    pol_angle,
    az_angle,
)
    s = 0.0
    k = length(coefs)
    l_max = Int(sqrt(k)) - 1
    spharm_array = real_spharm_array(l_max, pol_angle, az_angle)
    for j = 1:length(coefs)
        s += coefs[j] * spharm_array[j]
    end
    return s
end

function reconstruct_mesh(
    c_x::T,
    c_y::T,
    c_z::T,
    pol_angles,
    az_angles,
) where {T<:AbstractVector}
    @assert length(c_x) == length(c_y) == length(c_z)
    @assert length(pol_angles) == length(az_angles)
    n_points = length(pol_angles)

    x = zeros(Float64, n_points)
    y = zeros(Float64, n_points)
    z = zeros(Float64, n_points)
    for i = 1:n_points
        x[i] = reconstruct_from_spharm_coefs(c_x, pol_angles[i], az_angles[i])
        y[i] = reconstruct_from_spharm_coefs(c_y, pol_angles[i], az_angles[i])
        z[i] = reconstruct_from_spharm_coefs(c_z, pol_angles[i], az_angles[i])
    end
    return x, y, z
end

function reconstruct_mesh(
    c_x::T,
    c_y::T,
    c_z::T;
    div = 3,
) where {T<:AbstractVector}
    @assert length(c_x) == length(c_y) == length(c_z)
    sphere = icosphere(div, 1.0)
    verts = decompose(Point3f0, sphere)
    triangles = decompose(Face{3,Int}, sphere)
    verts_spherical = cart2sph.(verts)
    pol_angles = [v[1] for v in verts_spherical]
    az_angles = [v[2] for v in verts_spherical]
    x, y, z = reconstruct_mesh(c_x, c_y, c_z, pol_angles, az_angles)
    for i in eachindex(verts)
        verts[i] = Point3f0(x[i], y[i], z[i])
    end
    mesh = Mesh(verts, faces(sphere))
    return mesh
end

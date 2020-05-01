function rotation_matrix(theta, ax)
    x,y,z = LinearAlgebra.normalize(ax)
    c = cos(theta)
    s = sin(theta)
    R = [(c+x^2*(1-c)) (x*y*(1-c)-s*z) (z*x*(1-c)+s*y);
         (x*y*(1-c)+s*z) (c+y^2*(1-c)) (z*y*(1-c)-s*x);
         (x*z*(1-c)-s*y) (y*z*(1-c)+s*x) (c+z^2*(1-c))]
    return R
end
function parallel_transport(coordinate::Matrix{T_Number}) where {T_Number<:Number}
    r = copy(coordinate)
    n_points = size(r, 1)
    if size(r, 2) == 2
        r = cat(r, zeros(eltype(r), size(r, 1)), dims = 2)
    end
    t = 1.0:size(r,1)
    dt = (maximum(t) - minimum(t)) / length(r)
    dr = gradient(r)
    ddr = gradient(dr)
    dddr = gradient(ddr)


    T = similar(dr)
    N = similar(dr)
    B = similar(dr)
    T[1, :] = LinearAlgebra.normalize(dr[1, :])
    B[1, :] = LinearAlgebra.normalize([T[1,2], -T[1,1], 0 ])
    N[1, :] = LinearAlgebra.normalize(cross(B[1, :], T[1, :]))
    for i = 2:n_points
        T[i, :] = LinearAlgebra.normalize(dr[i, :])
        b_ = cross(T[i-1,:],T[i,:])
        if isapprox(norm(b_),0)
            B[i, :] = B[i-1,:]
        else
            b_ = LinearAlgebra.normalize(b_)
            phi = acos(dot(T[i-1,:],T[i,:]))
            R = rotation_matrix(phi,b_)
            B[i, :] = R * B[i-1,:]
        end
        N[i, :] = LinearAlgebra.normalize(cross(B[i, :], T[i, :]))
    end


    ds = Array{eltype(dr),1}(undef, n_points)
    curvature = similar(ds)
    torsion = similar(ds)

    for i in eachindex(ds)
        selected_dr = dr[i, :]
        selected_ddr = ddr[i, :]
        selected_dddr = dddr[i, :]
        ds[i] = norm(selected_dr)

        dr_cross_ddr = cross(selected_dr, selected_ddr)
        curvature_num = norm(dr_cross_ddr) # || dr X ddr ||
        curvature_denom = norm(selected_dr)^3 # || dr ||^3
        curvature[i] = curvature_num / curvature_denom

        torsion_num = dot(dr_cross_ddr, selected_dddr) # (dr X ddr) . dddr
        torsion_denom = norm(dr_cross_ddr)^2 # || dr X ddr ||^2
        torsion[i] = torsion_num / torsion_denom
    end
    speed = ds ./ dt
    return T,N,B,curvature,torsion,speed
end

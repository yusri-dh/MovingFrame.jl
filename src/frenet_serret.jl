function frenet_serret(coordinate::Matrix{T_Number}) where {T_Number<:Number}
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

    for i = 1:n_points

        T[i, :] = LinearAlgebra.normalize(dr[i, :])

        B[i, :] = LinearAlgebra.normalize(cross(dr[i, :], ddr[i, :]))

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

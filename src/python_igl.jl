const scipy_opt = PyNULL()



const igl = PyNULL()
const scipy = PyNULL()


# igl = pyimport("igl")
# scipy = pyimport("scipy")
function igl_gaussian_curvature(verts,faces)
    println("new gaussian curv")
    v = array_of_array_to_matrix(verts)
    f = array_of_array_to_matrix(faces) .-1
    K = igl.gaussian_curvature(v, f)
    m = igl.massmatrix(v, f, igl.MASSMATRIX_TYPE_BARYCENTRIC)
    # m = igl.massmatrix(v, f, igl.MASSMATRIX_TYPE_VORONOI)
    minv = scipy.sparse.diags(1.0 ./ m.diagonal())
    kn = minv.dot(K)
    return kn
end
function igl_mean_curvature(verts,faces)
    v = array_of_array_to_matrix(verts)
    f = array_of_array_to_matrix(faces) .-1
    l = igl.cotmatrix(v, f)
    m = igl.massmatrix(v, f, igl.MASSMATRIX_TYPE_BARYCENTRIC)
    # m = igl.massmatrix(v, f, igl.MASSMATRIX_TYPE_VORONOI)

    minv = scipy.sparse.diags(1.0 ./ m.diagonal())

    hn = -minv.dot(l.dot(v))
    h = scipy.linalg.norm(hn, axis=1)
    return h
end

function igl_mean_curvature_flow(v,f,n_step = 10, delta = 0.001)
    new_v = v .- centroid(v)
    new_v = new_v / sqrt(surface_area(new_v,f))
    l = igl.cotmatrix(new_v, f .- 1)


    for i in 1:n_step
        m = igl.massmatrix(new_v,f .- 1, igl.MASSMATRIX_TYPE_BARYCENTRIC)
        s = (m - delta * l)
        b = m.dot(new_v)
        new_v = scipy.sparse.linalg.spsolve(s, b)
        new_v = new_v .- centroid(new_v)
        new_v = new_v / sqrt(surface_area(new_v,f))

    end
    return new_v
end

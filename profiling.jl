using GeometryBasics;
using FileIO;
using LinearAlgebra;
using GaussianProcesses;
using MovingFrame;
import Makie;
using OwnTime;
using Profile;
const mv = MovingFrame
##11
meshfile = Base.download("https://github.com/libigl/libigl-tutorial-data/raw/master/bumpy-cube.obj")
shape = FileIO.load(FileIO.File{FileIO.format"OBJ"}(meshfile))
verts = decompose(Point3f0, shape)
faces = decompose(TriangleFace{Int}, shape)
vmat = MovingFrame.array_of_array_to_matrix(verts)
fmat = MovingFrame.array_of_array_to_matrix(faces)
##

function mean_curvature_flow_(verts, faces; n_step=10, delta=0.001,min_diff=1e-13)
    new_verts, _ = mv.volume_normalizing(mv.mesh_centering(verts, faces)...)
    v = MovingFrame.array_of_array_to_matrix(new_verts)
    f = MovingFrame.array_of_array_to_matrix(faces)
    new_v = v .- mv.centroid(v)
    new_v = new_v / sqrt(mv.surface_area2(new_v, f))
    l = mv.cotangent_laplacian(new_v, f)

    for i in 1:n_step
        old_v = copy(new_v)
        m = mv.mass_matrix_barycentric2(new_v, f)
        s = (m - delta * l)
        b = m * new_v
        new_v = s \ b
        new_v = new_v .- mv.centroid(new_v)
        new_v = new_v / sqrt(mv.surface_area2(new_v, f))
        diff = tr((new_v - old_v)' * mv.mass_matrix_barycentric2(new_v, f) * (new_v - old_v))
        # @show diff
        if diff <= min_diff
            println("$i,Converged")
            break
        end
    end
    new_verts::typeof(verts) = MovingFrame.matrix_to_vertex(new_v)
    return new_verts
end


function mean_curvature_flow_2(v::Matrix{Float32}, f::Matrix{Int64}; n_step=10, delta=0.001,min_diff=1e-13)
    new_v = v .- mv.centroid(v)
    new_v = new_v / sqrt(mv.surface_area2(new_v, f))
    l = mv.cotangent_laplacian(new_v, f)

    for i in 1:n_step
        old_v = copy(new_v)
        m = mv.mass_matrix_barycentric2(new_v, f)
        s = (m - delta * l)
        b = m * new_v
        new_v = s \ b
        new_v = new_v .- mv.centroid(new_v)
        new_v = new_v / sqrt(mv.surface_area2(new_v, f))
        diff = tr((new_v - old_v)' * mv.mass_matrix_barycentric2(new_v, f) * (new_v - old_v))
        # @show diff
        if diff <= min_diff
            println("Converged")
            break
        end
    end
    new_verts::typeof(verts) = MovingFrame.matrix_to_vertex(new_v)
    return new_verts
end
##
# spherization(shape)
Profile.clear()
# @profile mean_curvature_flow_(verts, faces,n_step=100)

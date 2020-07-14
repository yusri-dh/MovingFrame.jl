##
using GeometryBasics;
using FileIO;
using LinearAlgebra;
using GaussianProcesses;
using MovingFrame;
import Plots;
using Makie;
const plt = Plots
##
struct MovingCells
    coordinates::Matrix{Float64}
    time::Vector{Float64}
    smooth_n_points::Int64
    smooth_coordinates::Matrix{Float64}
    smooth_time::Vector{Float64}
    T::Matrix{Float64}
    N::Matrix{Float64}
    B::Matrix{Float64}
    curvature::Vector{Float64}
    torsion::Vector{Float64}
    speed::Vector{Float64}
    original_shapes::Vector{Mesh}
    reoriented_shapes::Vector{Mesh}
    Spharm_Cx::Matrix{Float64}
    Spharm_Cy::Matrix{Float64}
    Spharm_Cz::Matrix{Float64}
end

const base_dir = pwd();
const input_data_dir = joinpath(base_dir, "cell_shape_final/")
const input_files = readdir(input_data_dir);
const n_points = length(input_files)
##
coordinates = Matrix{Float64}(undef, n_points, 3)
time = collect(1.0:n_points)
shapes = Vector{Mesh}(undef, n_points)

for i = 1:length(input_files)
    shape = load(joinpath(input_data_dir, input_files[i]))
    shapes[i] = shape
    verts = decompose(Point3f0, shape)
    coordinates[i, :] .= reduce(+, verts) / length(verts)
end
shapes = (volume_normalizing ∘ mesh_centering).(shapes)
new_coordinates = coordinates .- coordinates[1, :]'
x, y, z = new_coordinates[:, 1], new_coordinates[:, 2], new_coordinates[:, 3]
##
gp_model_x = GP(time, x, MeanZero(), SE(4.0, 0.0))
gp_model_y = GP(time, y, MeanZero(), SE(4.0, 0.0))
gp_model_z = GP(time, z, MeanZero(), SE(4.0, 0.0))
smooth_x, _ = predict_y(gp_model_x, time)
smooth_y, _ = predict_y(gp_model_y, time)
smooth_z, _ = predict_y(gp_model_z, time)
smooth_coordinates = zeros(Float64, length(smooth_x), 3)
smooth_coordinates[:, 1] .= smooth_x
smooth_coordinates[:, 2] .= smooth_y
smooth_coordinates[:, 3] .= smooth_z

ori_index = round.(Int, time |> collect)
T, N, B, curvature, torsion, speed =
    parallel_transport(smooth_coordinates)
##
reoriented_shape = [
    mesh_reorientation(shapes[i], T[i, :], N[i, :], B[i, :])
    for i in eachindex(shapes)
]
spheres = MovingFrame.spherization.(reoriented_shape, 100)
##
spharm_coefs = [
    MovingFrame.spharm_descriptor(
        reoriented_shape[i],
        GeometryBasics.coordinates(spheres[i]),
    ) for i in eachindex(reoriented_shape)
]
Cx = MovingFrame.array_of_array_to_matrix([coef[:Cx] for coef in spharm_coefs])
Cy = MovingFrame.array_of_array_to_matrix([coef[:Cy] for coef in spharm_coefs])
Cz = MovingFrame.array_of_array_to_matrix([coef[:Cz] for coef in spharm_coefs])
##
xy_ecc = Cx[:, 4] ./ Cy[:, 2]
xz_ecc = Cx[:, 4] ./ Cz[:, 3]
yz_ecc = Cy[:, 2] ./ Cz[:, 3]
moving_cell = MovingCells(
    new_coordinates,
    time,
    length(time),
    smooth_coordinates,
    time,
    T,
    N,
    B,
    curvature,
    torsion,
    speed,
    shapes,
    reoriented_shape,
    Cx,
    Cy,
    Cz,
)
##
plt.plot(speed)
plt.plot(curvature)
plt.plot(torsion)

plt.plot([xy_ecc xz_ecc yz_ecc])

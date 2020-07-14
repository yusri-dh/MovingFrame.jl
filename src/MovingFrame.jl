module MovingFrame
using LinearAlgebra
using Statistics: mean, var, std
using GaussianProcesses
using GeometryBasics
using GSL
using PyCall
using SparseArrays
using DataFrames
using CSV
using FileIO: load, save

function __init__()
    copy!(igl, pyimport_conda("igl", "igl"))
    copy!(scipy, pyimport_conda("scipy", "scipy"))

end

include("./utils_movement.jl")
include("./smoothing.jl")
export Smoother, GaussianProcessSmoother, smoothing, fit, predict
include("./frenet_serret.jl")
export frenet_serret
include("./parallel_transport.jl")
export parallel_transport

include("spharm.jl")
export Plm,
    sph_Plm,
    Plm_array,
    Plm_array_index,
    real_spharm,
    real_spharm_array,
    j2lm,
    lm2j,
    spharm_coefs
include("shape_descriptor.jl")
export spharm_descriptor,
    reconstruct_from_spharm_coefs, reconstruct_mesh, sphere_parameterization
include("utils_shape.jl")
export cart2sph,
    sph2cart,
    mesh_centering,
    volume_normalizing,
    mesh_volume,
    create_mesh,
    centroid,
    internal_angles,
    double_surface_area,
    surface_area

include("sphere.jl")
include("python_igl.jl")
export igl_gaussian_curvature, igl_mean_curvature
include("curvature.jl")
export mass_matrix_barycentric,
    mass_matrix_barycentric2,
    cotangent_laplacian,
    gaussian_curvature,
    mean_curvature,
    spherization
include("moving_cells.jl")
export mesh_reorientation,
    MovingCells, get_cell, add_moving_frame_data!, add_shape_data!

include("fileIO.jl")
export load_mesh, save_mesh, create_group_df

function _save_shapes!(shapes::Vector{<:Mesh},dir_name)
    for i in eachindex(shapes)
        fname = "t" * string(i) * ".ply"
        shape=shapes[i]
        save_mesh(joinpath(dir_name, fname), shape)
    end
    0
end


function main(data_folder::String, track_id::Int, t_start::Int, t_end::Int)
    t0 = convert(Float64, t_start)
    t1 = convert(Float64, t_end)
    df = create_group_df(data_folder)
    cell = get_cell(df, track_id; t0 = t0, t1 = t1)
    coordinate = cell.coordinates
    length_step = 0.01
    add_moving_frame_data!(cell, step = length_step)
    if t0 == t1 == -1.0
        selected_df = filter(row -> row[:trackid] == track_id, df)
    else
        @assert t1 > t0 >= 1.0
        selected_df =
            filter(row -> (row.trackid == track_id) & (t0 <= row.t <= t1), df)
    end
    shape_folder = joinpath(data_folder,"shape")
    add_shape_data!(cell, selected_df, data_dir = shape_folder)
    str_t0::String = string(t0)
    str_t1::String = string(t1)
    str_id::String = string(track_id)

    dir_name = "movement_analysis" * "_" * str_id * "_" * str_t0 * "-" * str_t1
    dir_name = joinpath(data_folder,dir_name)
    if !isdir(dir_name)
        mkdir(dir_name)
    end
    movement_df = DataFrame(
        x = cell.smooth_coordinates[:, 1],
        y = cell.smooth_coordinates[:, 2],
        z = cell.smooth_coordinates[:, 3],
        time = cell.smooth_time,
        T = matrix_to_array_of_array(cell.T),
        N = matrix_to_array_of_array(cell.N),
        B = matrix_to_array_of_array(cell.B),
        curvature = cell.curvature,
        torsion = cell.torsion,
        speed = cell.speed,
    )
    fname = joinpath(dir_name, "result.csv")
    CSV.write(fname, movement_df)


    ori_dir_name =
        "original_shapes" * "_" * str_id * "_" * str_t0 * "-" * str_t1
    ori_dir_name = joinpath(data_folder,ori_dir_name)
    if !isdir(ori_dir_name)
        mkdir(ori_dir_name)
    end
    reori_dir_name =
        "reoriented_shapes" * "_" * str_id * "_" * str_t0 * "-" * str_t1
    reori_dir_name = joinpath(data_folder,reori_dir_name)
    if !isdir(reori_dir_name)
        mkdir(reori_dir_name)
    end
    _save_shapes!(cell.original_shapes,ori_dir_name)
    _save_shapes!(cell.reoriented_shapes,reori_dir_name)
    return 0
end

export main
end # module

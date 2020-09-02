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
    copy!(igl, pyimport_conda("igl", "igl","conda-forge"))
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
    surface_area,
    mesh_reorientation

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

include("fileIO.jl")
export load_mesh, save_mesh, create_group_df


end # module

### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ ff81125a-086d-11eb-2226-eddf1955c141
begin
	# Activate environment. You do not need to do this if you are not using
	# separate environment
	push!(LOAD_PATH, "$(homedir())/pCloudDrive/julia_project")
	using Pkg
	Pkg.activate("/home/yusri/pCloudDrive/julia_project/MovingFrame/")
end

# ╔═╡ 5df5c832-ee72-11ea-20bf-7b42b7e25571
begin
	using GeometryBasics
	using FileIO
	using LinearAlgebra;
	using GaussianProcesses;
	using PlutoUI
	using MovingFrame;
	import Plots;
	const mva = MovingFrame
	const plt = Plots
end

# ╔═╡ 67beff68-ee73-11ea-3b77-fbf331f8ce9c
# Create MovingCells Type
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
    original_shapes::Vector{GeometryBasics.Mesh}
    reoriented_shapes::Vector{GeometryBasics.Mesh}
    Spharm_Cx::Matrix{Float64}
    Spharm_Cy::Matrix{Float64}
    Spharm_Cz::Matrix{Float64}
end

# ╔═╡ 7139c988-ee73-11ea-1127-25bfb372000d
begin
	const base_dir = pwd();
	# You need to extract the cell_shape_final.zip first
	const input_data_dir = joinpath(base_dir, "cell_shape_final/")
	const input_files = readdir(input_data_dir);
	const n_points = length(input_files)
end

# ╔═╡ 9ebba1d6-ee73-11ea-1de2-d78152ef6b5e
begin
	# Calculating the cell trajectory
	coordinate = Matrix{Float64}(undef, n_points, 3)
	time = collect(1.0:n_points)
	shapes = Vector{GeometryBasics.Mesh}(undef, n_points)
	for i = 1:length(input_files)
		shape = load(joinpath(input_data_dir, input_files[i]))
		shapes[i] = shape
		verts = decompose(Point3f0, shape)
		# The trajectory is the mass centre of the cell
		coordinate[i, :] .= reduce(+, verts) / length(verts)
	end
	# Normalizing volume and translation
	shapes = (volume_normalizing ∘ mesh_centering).(shapes)
	new_coordinate = coordinate .- coordinate[1, :]'
end

# ╔═╡ 7d1cf63e-ee83-11ea-33e9-651efeecaf00
begin
	exp_lengthscale = 3.4 # lengthscale = e^3.4
	exp_variance= 0.0 #variance = e^0.0
	
	# Gaussian Processs using Zero mean function and 
	# Squared exponential kernel (note that hyperparameters are on the log scale)
	smoother =  GaussianProcessSmoother(MeanZero(), SE(exp_lengthscale, exp_variance))
	smooth_coordinates = smoothing(new_coordinate, smoother, step=1.)


	smooth_time = range(1.0, stop=time[end], length=size(smooth_coordinates, 1)) |> collect
	T, N, B, curvature, torsion, speed = parallel_transport(smooth_coordinates)
end

# ╔═╡ 5b0d2618-0878-11eb-0fe0-69b9357b91e9
# Reoriented the cell shape to the Moving Frame basis
reoriented_shape = [
    mesh_reorientation(shapes[i], T[i, :], N[i, :], B[i, :])
    for i in eachindex(shapes)
]

# ╔═╡ 29c8a4b8-087a-11eb-3780-4fa5a2c12d0e
# Spherical parameterization
spheres = MovingFrame.spherization.(reoriented_shape, n_step=100)


# ╔═╡ ed0b2ee6-090d-11eb-3516-9b7788a1859b
begin
	# Calculating spherical harmonics coefficients
	spharm_coefs = [
	    MovingFrame.spharm_descriptor(
	        reoriented_shape[i],
	        GeometryBasics.coordinates(spheres[i]),
	    ) for i in eachindex(reoriented_shape)
	]
	Cx = MovingFrame.array_of_array_to_matrix([coef[:Cx] for coef in spharm_coefs])
	Cy = MovingFrame.array_of_array_to_matrix([coef[:Cy] for coef in spharm_coefs])
	Cz = MovingFrame.array_of_array_to_matrix([coef[:Cz] for coef in spharm_coefs])
end

# ╔═╡ 18462938-090e-11eb-2154-45881e7ef9f1
begin
	# Calculating shape eccentricity
	xy_ecc = Cx[:, 4] ./ Cy[:, 2]
	xz_ecc = Cx[:, 4] ./ Cz[:, 3]
	yz_ecc = Cy[:, 2] ./ Cz[:, 3]
end

# ╔═╡ cd0cd780-ee89-11ea-0e25-45fdbb85346d
begin
plt.plot(smooth_coordinates[:,1],smooth_coordinates[:,2],smooth_coordinates[:,3],label ="SMOOTH_TRAJECTORY")
plt.scatter!(new_coordinate[:,1],new_coordinate[:,2],new_coordinate[:,3],label ="OBSERVATION",markersize=1.0)

end

# ╔═╡ 46ff40b6-090e-11eb-1216-b1a6a3af5695
plt.plot(smooth_time,[xy_ecc xz_ecc],label=["XY Eccentricity" "XZ Eccentricity"])

# ╔═╡ 04cdf2de-3dcd-11eb-2123-9fef1cb5553e
begin
	p1 = plt.plot(smooth_time,speed,ylabel = "Speed")
	p2 = plt.plot(smooth_time,curvature,ylabel = "Curvature",)
	p3 = plt.plot(smooth_time,torsion,ylabel="Torsion",)
	plt.plot(p1,p2,p3,legend=false)
end

# ╔═╡ Cell order:
# ╠═ff81125a-086d-11eb-2226-eddf1955c141
# ╠═5df5c832-ee72-11ea-20bf-7b42b7e25571
# ╠═67beff68-ee73-11ea-3b77-fbf331f8ce9c
# ╠═7139c988-ee73-11ea-1127-25bfb372000d
# ╠═9ebba1d6-ee73-11ea-1de2-d78152ef6b5e
# ╠═7d1cf63e-ee83-11ea-33e9-651efeecaf00
# ╠═5b0d2618-0878-11eb-0fe0-69b9357b91e9
# ╠═29c8a4b8-087a-11eb-3780-4fa5a2c12d0e
# ╠═ed0b2ee6-090d-11eb-3516-9b7788a1859b
# ╠═18462938-090e-11eb-2154-45881e7ef9f1
# ╠═cd0cd780-ee89-11ea-0e25-45fdbb85346d
# ╠═46ff40b6-090e-11eb-1216-b1a6a3af5695
# ╠═04cdf2de-3dcd-11eb-2123-9fef1cb5553e

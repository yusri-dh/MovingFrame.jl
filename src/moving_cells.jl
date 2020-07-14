mutable struct MovingCells
    track_id::Int64
    coordinates::Array{Float64,2}
    time::Array{Float64,1}
    smooth_n_points::Int64
    smooth_coordinates::Array{Float64,2}
    smooth_time::Array{Float64,1}
    T::Array{Float64,2}
    N::Array{Float64,2}
    B::Array{Float64,2}
    curvature::Array{Float64,1}
    torsion::Array{Float64,1}
    speed::Array{Float64,1}
    original_shapes::Array{Mesh,1}
    reoriented_shapes::Array{Mesh,1}
    Spharm_Cx::Array{Float64,2}
    Spharm_Cy::Array{Float64,2}
    Spharm_Cz::Array{Float64,2}

    function MovingCells()
        return new()
    end

    function MovingCells(
        track_id::Int64,
        coordinates::Array{Float64,2},
        time::Array{Float64,1},
    )
        cell = new()
        cell.track_id = track_id
        cell.coordinates = coordinates .- coordinates[1, :]'
        cell.time = time

        return cell
    end
end



function get_cell(
    df::DataFrame,
    track_id::Int;
    t0 = -1.0,
    t1 = -1.0,
)::MovingCells

    if t0 == t1 == -1.0
        selected_df = filter(row -> row[:trackid] == track_id, df)
    else
        @assert t1 > t0 >= 1.0
        selected_df =
            filter(row -> (row.trackid == track_id) & (t0 <= row.t <= t1), df)
    end
    coordinates = convert(Matrix, selected_df[!, [:x, :y, :z]])
    time = collect(1.0:size(coordinates, 1))
    cell = MovingCells(track_id, coordinates, time)
    return cell
end




function add_moving_frame_data!(cell::MovingCells; step = 0.01)
    coordinates = cell.coordinates
    smooth_coordinates = smoothing(coordinates, step = step)
    cell.smooth_coordinates = smooth_coordinates
    cell.smooth_n_points = size(smooth_coordinates, 1)
    cell.smooth_time =
        range(minimum(cell.time), stop = maximum(cell.time), step = step)
    T, N, B, curvature, torsion, speed =
        parallel_transport(smooth_coordinates)
    cell.T, cell.N, cell.B = T, N, B
    cell.curvature, cell.torsion, cell.speed = curvature, torsion, speed
end


function extract_shape_data(df; data_dir = pwd())
    @assert length(unique(df.trackid)) == 1
    cell_id = df.cellid
    gc_cell_id = df.gc_cellid
    t = convert.(Int, df.t)
    shapes = Vector{Mesh}(undef, length(gc_cell_id))
    for i in eachindex(gc_cell_id)
        fname = joinpath(
            data_dir,
            "t" * lpad(t[i], 4, "0"),
            "obj0",
            "cell" * lpad(gc_cell_id[i], 5, "0") * ".obj",
        )
        # @show fname
        shape = load(fname) |> mesh_centering |> volume_normalizing

        shapes[i] = load(fname)
    end
    return shapes
end

function add_shape_data!(cell::MovingCells,df;data_dir = pwd())
    shapes = extract_shape_data(df,data_dir=data_dir)
    cell.original_shapes = shapes
    smooth_time = cell.smooth_time
    step = smooth_time[2] - smooth_time[1]
    selected_time = [
                    1 + round(Int, 1 / step) * x
                    for
                    x in range(
                        0,
                        stop = div(length(smooth_time), round(Int, 1 / step)),
                    )
                ]
    T = cell.T[selected_time, :]
    N = cell.N[selected_time, :]
    B = cell.B[selected_time, :]
    reoriented_shapes = [mesh_reorientation(shapes[i], T[i, :], N[i, :], B[i, :]) for i in eachindex(shapes)]
    cell.reoriented_shapes = reoriented_shapes
end



function add_spharm_coef_data!(cell::MovingCells)
    reoriented_shape = cell.reoriented_shapes
    spheres = sh.spherization.(reoriented_shape, 100)
    spharm_coefs = [
        sh.spharm_descriptor(reoriented_shape[i], spheres[i].vertices)
        for i in eachindex(reoriented_shape)
    ]
    cell.Spharm_Cx =
        sh.array_of_array_to_matrix([coef[:Cx] for coef in spharm_coefs])
    cell.Spharm_Cy =
        sh.array_of_array_to_matrix([coef[:Cy] for coef in spharm_coefs])
    cell.Spharm_Cz =
        sh.array_of_array_to_matrix([coef[:Cz] for coef in spharm_coefs])
end

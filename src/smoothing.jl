# include("./utils.jl")

abstract type Smoother end
mutable struct GaussianProcessSmoother <: Smoother
    smoother_function :: Function
    kernel
    mean_function
    GaussianProcessSmoother(x,y,z) = GaussianProcessSmoother(x,y,z)
    GaussianProcessSmoother() = new()
    function GaussianProcessSmoother(mean_function,kernel)
        gp = new()
        gp.smoother_function = GP
        gp.mean_function = mean_function
        gp.kernel = kernel
        return gp
    end
end

function fit(smoother::GaussianProcessSmoother,x,y)
    model = smoother.smoother_function(
                        x,
                        y,
                        smoother.mean_function,
                        smoother.kernel
                        )
    return model
end

function predict(model::GPBase,smoothed_t)
    smoothed_x,_ = predict_y(model, smoothed_t)
    return smoothed_x
end

function smoothing(x::AbstractVector,smoother::Smoother;step=0.01)
    t = 1.0:length(x)
    standardized_x = standardize(x)
    model = fit(smoother,t,standardized_x)
    lower_bound, upper_bound = minimum(t), maximum(t)
    smoothed_t = range(lower_bound, stop = upper_bound, step = step)
    n_points = length(smoothed_t)
    smoothed_x = predict(model, smoothed_t)
    smoothed_x = reverse_standardization(smoothed_x,x)
    return smoothed_x
end
function smoothing(x::AbstractVector;step=0.01)
    smoother =  GaussianProcessSmoother(MeanZero(),SE(0.0,0.0))
    smoothed_x = smoothing(x,smoother,step=step)
    return smoothed_x
end

function smoothing(x::Matrix{T},smoothers::Vector{Smoother};step=0.01) where {T<:Number}
    n_row, n_col = size(x,1),size(x,2)
    @assert n_col == length(smoothers)
    t = 1.0:n_row
    lower_bound, upper_bound = minimum(t), maximum(t)
    smoothed_t = range(lower_bound, stop = upper_bound, step = step)
    n_points = length(smoothed_t)
    new_x = zeros(eltype(x),n_points,n_col)

    for i in 1:n_col
        smoothed_x = smoothing(x[:,i],smoothers[i],step=step)
        new_x[:, i] = smoothed_x
    end
    return new_x
end
function smoothing(x::Matrix{T},smoother::Smoother;step=0.01) where {T<:Number}
    smoother_x = smoother
    smoother_y = deepcopy(smoother)
    smoother_z = deepcopy(smoother)
    smoothers = [smoother_x,smoother_y,smoother_z]
    return smoothing(x,smoothers;step=step)
end
function smoothing(x::Matrix{T};step=0.01) where {T<:Number}
    n_row, n_col = size(x,1),size(x,2)
    t = 1.0:n_row
    lower_bound, upper_bound = minimum(t), maximum(t)
    smoothed_t = range(lower_bound, stop = upper_bound, step = step)
    n_points = length(smoothed_t)
    new_x = zeros(eltype(x),n_points,n_col)
    for i in 1:n_col
        smoothed_x = smoothing(x[:,i],step=step)
        new_x[:, i] = smoothed_x
    end
    return new_x
end

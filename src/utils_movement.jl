
function standardize(x::AbstractVector)
    std_x = std(x)
    if std_x == zero(std_x)
        return zero(x)
    else
        return (x .- mean(x)) ./ std_x
    end
end

function standardize(x::AbstractMatrix;dims = 1 )
    standardized_x = similar(x)
    if dims == 1
        for i in 1:size(x,2)
            standardized_x[:,i] .= standardize(x[:,i])
        end
    elseif dims == 2
        for i in 1:size(x,1)
            standardized_x[i,:] .= standardize(x[i,:])
        end
    end
    return standardized_x
end
function reverse_standardization(standardized_x::T, x::T) where T<: AbstractVector
    std_x =std(x)
    if std_x == zero(std_x)
        temp = x[1]
        return fill(temp,size(standardized_x))
    else
        return (standardized_x .* std(x)) .+ mean(x)
    end
end

function reverse_standardization(standardized_x::T, x::T; dims=1) where T<: AbstractMatrix
    new_x =similar(standardized_x)
    if dims == 1
        for i in 1:size(x,2)
            new_x[:,i] .= reverse_standardization(standardized_x[:,i],x[:,i])
        end
    elseif dims == 2
        for i in 1:size(x,1)
            new_x[i,:] .= reverse_standardization(standardized_x[i,:],x[i,:])
        end
    end
    return new_x
end

function gradient(x::Vector{T}) where {T<:Number}
    res = similar(x)
    for i = 1:length(x)
        if i == 1
            res[i] = x[i+1] - x[i]
        elseif i == size(x, 1)
            res[i] = x[i] - x[i-1]
        else
            res[i] = (x[i+1] - x[i-1]) / 2
        end
    end
    return res
end

function gradient(m::Matrix{T}) where {T<:Number}
    res = similar(m)
    for i = 1:size(m, 1)
        if i == 1
            res[i, :] = m[i+1, :] - m[i, :]
        elseif i == size(m, 1)
            res[i, :] = m[i, :] - m[i-1, :]
        else
            res[i, :] = (m[i+1, :] - m[i-1, :]) / 2
        end
    end
    return res
end


# function TNB_reorientation(verts::AbstractVector, T, N, B)
#     new_verts = similar(verts)
#     for i in eachindex(new_verts)
#         new_verts[i] = [T N B] \ verts[i]
#     end
#     return new_verts
# end

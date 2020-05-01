# using LinearAlgebra
# include("sphere.jl")
##

"""
This function compute the associated Legendre polynomial P(l,m,x)
    for m >= 0, l >= m, and abs(x) <= 1
"""
function Plm(l::Int, m::Int, x::AbstractFloat)
    @assert 0 <= m <= l
    pmm = 1.0
    if m > 0
        somx2 = sqrt((1.0 - x) * (1.0 + x))
        fact = 1.0
        for i = 1:m
            pmm *= (-fact) * somx2
            fact += 2.0
        end
    end

    if l == m
        return pmm
    end
    pmmp1 = x * (2.0 * m + 1.0) * pmm
    if (l == (m + 1))
        return pmmp1
    end

    pll = 0.0

    for ll = (m+2):l
        pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
    end

    return pll
end
"""
This function compute the normalized associated Legendre polynomial P(l,m,x)
    for m >= 0, l >= m, and abs(x) <= 1. It is alias  of GSL sf_legendre_Plm
"""
function sph_Plm(l::Int, m::Int, x::AbstractFloat)
    @assert 0 <= m <= l
    return sf_legendre_sphPlm(l, m, x)
end
"""
This function compute all of the normalized associated Legendre polynomial P(l,m)
    for m >= 0, l >= m, and abs(x) <= 1. It is alias  of GSL sf_legendre_array.
    This implementation include the Condon-Shortley phase factor.
"""
function Plm_array(l_max::Int, x::AbstractFloat; normal = :SPHARM)
    normalization = Dict(
        :SPHARM => GSL_SF_LEGENDRE_SPHARM,
        :NONE => GSL_SF_LEGENDRE_NONE,
        :SCHMIDT => GSL_SF_LEGENDRE_SCHMIDT,
        :FULL => GSL_SF_LEGENDRE_FULL,
    )
    return sf_legendre_array_e(normalization[normal], l_max, x,1)
end
"""
This function returns the index of Plm_array corresponding to P(l,m)
"""
function Plm_array_index(l::Int, m::Int)
    @assert 0 <= m <= l
    return sf_legendre_array_index(l,m) + 1
end
"""
This function returns the the general real spherical harmonics function Y(l,m) at
    pol angle and az_angle. l ∈ (0,..,n) and ∈ (-l,..l) m are the
    order and degree of spherical harmonics respectively.
    Pol_angle ∈ (0,..pi) and az_angle ∈ (0,2pi) are the polar angle and azimuthal angle
    coordinate.
"""
function real_spharm(l::Int, m::Int, pol_angle, az_angle)
    sqrt2 = sqrt(2.0)
    if m == 0
        return sph_Plm(l, m, cos(pol_angle))
    elseif m > 0
        return sqrt2 * cos(m * az_angle) * sph_Plm(l, m, cos(pol_angle))
    elseif m < 0
        return sqrt2 * sin(-m * az_angle) * sph_Plm(l, -m, cos(pol_angle))
    end
end
"""
This function returns all of general real spherical harmonics function Y(l,m)
    with order l.
    Pol_angle ∈ (0,..pi) and az_angle ∈ (0,2pi) are the polar angle and azimuthal angle
    coordinate.
"""
function real_spharm_array(l_max::Int, pol_angle, az_angle)
    sqrt2 = sqrt(2.0)
    len_array = (l_max + 1)^2
    plm_array = Plm_array(l_max, cos(pol_angle); normal = :SPHARM)
    spharm_array = Vector{Float64}(undef, len_array)
    for j = 1:len_array
        l, m = j2lm(j)
        idx = Plm_array_index(l, abs(m))
        if m == 0
            spharm_array[j] = plm_array[idx]
        elseif m > 0
            spharm_array[j] = sqrt2 * cos(m * az_angle) * plm_array[idx]
        elseif m < 0
            spharm_array[j] = sqrt2 * sin(-m * az_angle) * plm_array[idx]
        end
    end
    return spharm_array
end


"""
This function returns the order l and degree m of Y(l,m) from linear
    index j of spharm_array
"""
function j2lm(j)
    j = j - 1
    l = floor(Int, sqrt(j))
    m = j - l * (l + 1)
    return l, m
end
"""
This function returns the linear index j of spharm_array
    from the order l and degree m of Y(l,m)
"""
function lm2j(l, m)
    return l^2 + l + m + 1
end

"""
This function returns the The matrix of general real spherical
    harmonics function Y(l,m). The row index i  corresponds to points x_i and
    column index j obtained from j =  l^2 + l + m +1 where l and m are the order
    and degree of spherical harmonics respectively.
"""
function spharm_mat(l_max, pol_angle::T, az_angle::T) where {T<:AbstractVector}
    n_row = length(pol_angle)
    n_column = (l_max + 1)^2
    mat = Matrix{eltype(pol_angle)}(undef, n_row, n_column)
    for i = 1:n_row
        mat[i, :] .= real_spharm_array(l_max, pol_angle[i], az_angle[i])
    end
    return mat
end

"""
This function returns the The matrix of coefficient of SH function Y(l,m). The
    row index i  corresponds to points x_i and column index j obtained from j =  l^2 + l + m +1 where l and m are the order
    and degree of spherical harmonics respectively.
"""
function spharm_coefs(l_max::Int, x, pol_angle,az_angle)
    mat = spharm_mat(l_max, pol_angle, az_angle)
    coefs = mat \ x
    return coefs
end

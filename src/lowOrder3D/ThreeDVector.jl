#===============================================================================
    ThreeDVector

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDVector
    x :: Float64
    y :: Float64
    z :: Float64

    function ThreeDVector(x_ :: T1, y_ :: T2, z_ :: T3) where
        {T1 <: Real, T2 <: Real,  T3 <: Real}
        new(x_, y_, z_)
    end

    function ThreeDVector(array :: Vector{T}) where T <: Real
        @assert(size(array)[1] == 3)
        x = Float64(array[1])
        y = Float64(array[2])
        z = Float64(array[3])
        new(x, y, z)
    end
end

""" Convert type ThreeDVector to Array{Float64}"""
function convert(::Type{Array{T, 1}}, a::ThreeDVector) where T <: Real
    return [a.x, a.y, a.z]
end

""" Convert type Array{Real} to ThreeDVector"""
function convert(::Type{ThreeDVector}, a::Array{T, 1}) where T<: Real
    @assert(size(a)[1] == 3)
    b = ThreeDVector(a[1], a[2], a[3])
    return b
end

function Base.:+(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector([
        a.x + b.x,
        a.y + b.y,
        a.z + b.z ])
    return c
end

function Base.:-(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector([a.x - b.x, a.y - b.y, a.z - b.z])
    return c
end

function Base.:-(a::ThreeDVector)
    c = ThreeDVector(-a.x, -a.y , -a.z)
    return c
end

function Base.:*(a::ThreeDVector, b::T) where T <: Real
    c = ThreeDVector([
        a.x * b,
        a.y * b,
        a.z * b ])
    return c
end

function Base.:*(a::T, b::ThreeDVector) where T <: Real
    return b * a
end

function Base.:*(a::Matrix{Real}, b::ThreeDVector)
    @assert(size(matrix) == (3,3))
    c = ThreeDVector([
        dot(a[1,:], b),
        dot(a[2,:], b),
        dot(a[3,:], b)
    ])
    return c
end

function Base.:/(a::ThreeDVector, b::T) where T <: Real
    c = ThreeDVector([
        a.x / b,
        a.y / b,
        a.z / b ])
    return c
end

"""Cross product of two ThreeDVectors"""
function cross(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector(
        [a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x ])
    return c
end

"""Dot product of two ThreeDVectors"""
function dot(a::ThreeDVector, b::ThreeDVector)
    return a.x * b.x + a.y * b.y + a.z * b.z
end


function dot(a::Vector{Real}, b::ThreeDVector)
    @assert(length(a) == 3)
    return a[1] * b.x + a[2] * b.y + a[3] * b.z
end

function dot(a::ThreeDVector, b::Vector{Real})
    return dot(b, a)
end

function Base.abs(a::ThreeDVector)
    return sqrt(a.x^2 + a.y^2 + a.z^2)
end

""" Set vector a to zero"""
function zero!(a::ThreeDVector)
    a.x = 0.0
    a.y = 0.0
    a.z = 0.0
    return Void
end

"""Return vector of length 1 with same direction"""
function unit(a::ThreeDVector)
    b = abs(a)
    return a / b
end

"""Set a vector normalised to length 1 with same direction"""
function unit!(a::ThreeDVector)
    b = abs(a)
    a.x /= b
    a.y /= b
    a.z /= b
    return Void
end

function iszero(a::ThreeDVector)
    if a.x == 0.0 && a.y == 0.0 && a.z == 0.0
        return true
    else
        return false
    end
end

function Base.getindex(a::ThreeDVector, i::Int)  where Int <: Integer
    if i == 1 return a.x
    elseif i == 2 return a.y
    elseif i == 3 return a.z
    else throw(Core.BoundsError)
    end
end

function Base.size(a::ThreeDVector)
    return [3]
end

function Base.isfinite(a::ThreeDVector)
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z)
end

function rotate_about_x(point::ThreeDVector, angle :: Float64)
    s = sin(angle)
    c = cos(angle)
    mtrw = [1 0 0; 0 c -s; 0 s c]
    return convert(ThreeDVector, mtrw * convert(Vector{Float64}, point))
end

function rotate_about_y(point::ThreeDVector, angle :: Float64)
    s = sin(angle)
    c = cos(angle)
    mtrw = [c 0 s; 0 1 0; -s 0 c]
    return convert(ThreeDVector, mtrw * convert(Vector{Float64}, point))
end

function rotate_about_z(point::ThreeDVector, angle :: Float64)
    s = sin(angle)
    c = cos(angle)
    mtrw = [c -s 0; s c 0; 0 0 1]
    return convert(ThreeDVector, mtrw * convert(Vector{Float64}, point))
end

function rotate_about_axis(
    point::ThreeDVector,
    axis :: ThreeDVector,
    angle :: Float64)
    a = unit(axis)
    e = a.x
    f = a.y
    g = a.z
    s = sin(angle)
    c = cos(angle)
    mtrw = [e*e*(1-c)+c     e*f*(1-c)-g*s   g*e*(1-c)+f*s;
            e*f*(1-c)+f*c   f*f*(1-c)+c     g*f*(1-c)-e*s;
            e*g*(1-c)-g*c   g*f*(1-c)+e*s   g*g*(1-c)+c     ]
    return convert(ThreeDVector, mtrw * convert(Vector{Float64}, point))
end

#= END ThreeDVector ----------------------------------------------------------=#

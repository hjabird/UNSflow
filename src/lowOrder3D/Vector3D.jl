#===============================================================================
    Vector3D.jl

    Representation of a three dimensional vector.

    Why?:   Julia does not implement statically sized vectors by default. Odd.

    Initial code: HJAB 2018

    Copyright (c) 2018 HJA Bird

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to
    deal in the Software without restriction, including without limitation the
    rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
    sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
    IN THE SOFTWARE.
------------------------------------------------------------------------------=#
import ForwardDiff

mutable struct Vector3D
    x :: Float64
    y :: Float64
    z :: Float64

    function Vector3D(x_ :: T1, y_ :: T2, z_ :: T3) where
        {T1 <: Real, T2 <: Real,  T3 <: Real}
        new(Float64(x_), Float64(y_), Float64(z_))
    end

    function Vector3D(array :: Vector{T}) where T <: Real
        @assert(size(array)[1] == 3)
        x = Float64(array[1])
        y = Float64(array[2])
        z = Float64(array[3])
        new(x, y, z)
    end

    function Vector3D(a::Vector3D)
        new(a.x, a.y, a.z)
    end
end

function Base.convert(::Type{Vector{T}}, a::Vector3D) where T <: Real
    return [a.x, a.y, a.z]
end

function Base.convert(::Type{Vector3D}, a::Vector{T}) where T<: Real
    @assert(size(a)[1] == 3)
    b = Vector3D(a[1], a[2], a[3])
    return b
end

function Base.convert(::Type{Vector{T}}, a::Vector{Vector3D}) where T <: Real
    arr = Array{T, 1}(undef, 3*length(a))
    for i = 1:length(a)
        arr[i*3 - 2] = a[i].x
        arr[i*3 - 1] = a[i].y
        arr[i*3 - 0] = a[i].z
    end
    return arr
end

function Base.convert(::Type{Vector{Vector3D}}, a::Array{T, 1}) where T <: Real
    @assert(length(a) % 3 == 0, string("Conversion of Array{T<:Real,1} to"*
        " Array{UNSflow.Vector3D} requires that the array's length be a "*
        " multiple of 3. Length was ", length(a)))
    arr = Vector{Vector3D}(undef, length(a)/3)
    for i = 1:length(arr)
        arr[i].x = arr[i*3 - 2]
        arr[i].y = arr[i*3 - 1]
        arr[i].z = arr[i*3 - 0]
    end
    return arr
end

function Base.convert(::Type{Matrix{T}}, a::Vector{Vector3D}) where T <: Real
    arr = Matrix{T}(undef, (3, length(a)))
    for i = 1:length(a)
        arr[:, i] = [a[i].x, a[i].y, a[i].z]
    end
    return arr
end

function Base.convert(::Type{Vector{Vector3D}}, a::Matrix{T}) where T <: Real
    @assert(size(a)[1] == 3, string("Expected size(a::Matrix)[1] to be equal "*
        " 3. size(a::Matrix) = ", size(a), size(a)[2] == 3 ? ". The second "*
        "dimension == 3. To transpose the matrix, use the apostrophe (')"*
        " operator." : "."))
    arr = Vector{Vector3D}(undef, size(a)[2])
    for i = 1:length(arr)
        arr[i] = convert(Vector3D, a[:, i])
    end
    return arr
end

function Base.:+(a::Vector3D, b::Vector3D)
    c = Vector3D([
        a.x + b.x,
        a.y + b.y,
        a.z + b.z ])
    return c
end

function Base.:-(a::Vector3D, b::Vector3D)
    c = Vector3D([a.x - b.x, a.y - b.y, a.z - b.z])
    return c
end

function Base.:-(a::Vector3D)
    c = Vector3D(-a.x, -a.y , -a.z)
    return c
end

function Base.:*(a::Vector3D, b::T) where T <: Real
    c = Vector3D([
        a.x * b,
        a.y * b,
        a.z * b ])
    return c
end

function Base.:*(a::T, b::Vector3D) where T <: Real
    return b * a
end

function Base.:*(a::Array{T, 2}, b::Vector3D) where {T <: Real}
    @assert(size(a) == (3,3))
    c = Vector3D([
        dot(a[1,:], b),
        dot(a[2,:], b),
        dot(a[3,:], b)
    ])
    return c
end

function Base.:/(a::Vector3D, b::T) where T <: Real
    c = Vector3D([
        a.x / b,
        a.y / b,
        a.z / b ])
    return c
end

function cross(a::Vector3D, b::Vector3D)
    c = Vector3D(
        [a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x ])
    return c
end

function dot(a::Vector3D, b::Vector3D)
    return a.x * b.x + a.y * b.y + a.z * b.z
end


function dot(a::Array{T, 1}, b::Vector3D) where T <: Real
    @assert(length(a) == 3)
    return a[1] * b.x + a[2] * b.y + a[3] * b.z
end

function dot(a::Vector3D, b::Vector{Real})
    return dot(b, a)
end

function Base.abs(a::Vector3D)
    return sqrt(a.x^2 + a.y^2 + a.z^2)
end

function zero!(a::Vector3D)
    a.x = 0.0
    a.y = 0.0
    a.z = 0.0
    return Void
end

function unit(a::Vector3D)
    b = abs(a)
    return a / b
end

function unit!(a::Vector3D)
    b = abs(a)
    a.x /= b
    a.y /= b
    a.z /= b
    return
end

function Base.:(==)(a::Vector3D, b::Vector3D)
    if a.x == b.x && a.y == b.y && a.z == b.z
        return true
    else
        return false
    end
end

function iszero(a::Vector3D)
    if a.x == 0.0 && a.y == 0.0 && a.z == 0.0
        return true
    else
        return false
    end
end

function Base.getindex(a::Vector3D, i::Int)  where Int <: Integer
    if i == 1 return a.x
    elseif i == 2 return a.y
    elseif i == 3 return a.z
    else throw(Core.BoundsError)
    end
end

function Base.size(a::Vector3D)
    return [3]
end

function Base.length(a::Vector3D)
    return 3
end

function Base.iterate(a::Vector3D, state=(a, 1))
    element, count = state
    if count > 3
        return nothing
    end
    return (a[count], (a, count + 1))
end

function Base.eltype(::Type{Tuple{Vector3D, T}}) where T <: Int
    return Float64
end

function outer(a::Vector3D, b::Vector3D)
    return [
        a.x*b.x a.x*b.y a.x*b.z;
        a.y*b.x a.y*b.y a.y*b.z;
        a.z*b.x a.z*b.y a.z*b.z;
        ]
end

function Base.isfinite(a::Vector3D)
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z)
end

function rotate_about_x(point::Vector3D, angle :: Float64)
    s = sin(angle)
    c = cos(angle)
    mtrw = [1 0 0; 0 c -s; 0 s c]
    return convert(Vector3D, mtrw * convert(Vector{Float64}, point))
end

function rotate_about_y(point::Vector3D, angle :: Float64)
    s = sin(angle)
    c = cos(angle)
    mtrw = [c 0 s; 0 1 0; -s 0 c]
    return convert(Vector3D, mtrw * convert(Vector{Float64}, point))
end

function rotate_about_z(point::Vector3D, angle :: Float64)
    s = sin(angle)
    c = cos(angle)
    mtrw = [c -s 0; s c 0; 0 0 1]
    return convert(Vector3D, mtrw * convert(Vector{Float64}, point))
end

function rotate_about_axis(
    point::Vector3D,
    axis :: Vector3D,
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
    return convert(Vector3D, mtrw * convert(Vector{Float64}, point))
end
#= END Vector3D ----------------------------------------------------------=#

#===============================================================================
    BilinearQuad.jl

    Represents a bilinear quad element with local coordinates [-1,1], [-1,1]

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

mutable struct BilinearQuad <: DiscreteGeometry3D
    c1 :: Vector3D
    c2 :: Vector3D
    c3 :: Vector3D
    c4 :: Vector3D

    function BilinearQuad(a::Vector{Vector3D})
        @assert(length(a) == 4, "BilinearQuad is defined by four points")
        new(a[1], a[2], a[3], a[4])
    end
    function BilinearQuad(
        corner1::Vector3D, corner2::Vector3D,
        corner3::Vector3D, corner4::Vector3D)
        new(corner1, corner2, corner3, corner4)
    end
end

# Return map a local coordinate to a point in space
function evaluate(a::BilinearQuad, local_coord::Vector{T}) where T <: Real
    @assert(length(position) == 2, "BilinearQuad is two dimensionsal.")
    x = local_coord[1]
    y = local_coord[2]
    t1 = a.c1 * (x-1)(y-1) / 4
    t2 = - a.c2 * (x+1)(y-1) / 4
    t3 = a.c3 * (x+1)(y+1) / 4
    t4 = - a.c4 * (x-1)(y+1) / 4
    return t1 + t2 + t3 + t4
end

# Return a vector of coordinates the interpolation points of the geometry
function coords(a::BilinearQuad)
    return [a.c1, a.c2, a.c3, a.c4]
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::BilinearQuad)
    return 2
end

# Get the number of control points of an object
function number_of_control_points(a::BilinearQuad)
    return 4
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::BilinearQuad, position::Vector{T}) where T <: Real
    @assert(length(position) == 2, "BilinearQuad is two dimensionsal.")
    return (-1 <= position[1] <= 1) && (-1 <= position[2] <= 1)
end

"Make a bilinear quad planar using the normal at it centre."
function flatten!(a::BilinearQuad)

    return
end

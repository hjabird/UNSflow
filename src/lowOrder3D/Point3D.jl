#===============================================================================
    DiscreteGeometry3D.jl

    Represents a geometry defined by points in 3D cartesian space.

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

mutable struct Point3D <: DiscreteGeometry3D
    coord :: Vector3D

    function Point3D(a::Vector3D)
        new(a)
    end
end


# Return map a local coordinate to a point in space
function evaluate(a::Point3D, local_coord::Vector{Real})
    @assert(length(position) == 0, "Point3D is zero dimensionsal.")
    return a.coord
end

# Return a vector of coordinates the interpolation points of the geometry
function coords(a::Point3D)
    return [a.coord]
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::Point3D)
    return 0
end

# Get the number of control points of an object
function number_of_control_points(a::Point3D)
    return 1
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::Point3D, position::Vector{Real})
    @assert(length(position) == 0, "Point3D is zero dimensionsal.")
    return true
end

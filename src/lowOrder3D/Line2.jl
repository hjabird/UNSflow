#===============================================================================
    Line2.jl

    Represents a linear line in cartesian space

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

mutable struct Line2 <: DiscreteGeometry3D
    start_coord :: Vector3D
    end_coord :: Vector3D

    function Line2(a::Vector{Vector3D})
        @assert(length(a) == 2, "Line2 is defined by two points")
        new(a[1], a[2])
    end
    function Line2(start::Vector3D, ending::Vector3D)
        new(start, ending)
    end
end

# Return map a local coordinate to a point in space
function evaluate(a::Line2, local_coord::Vector{Real})
    @assert(length(position) == 1, "Line2 is one dimensionsal.")
    return 0.5 * (a.start_coord * (local_coord[1] - 1) +
        a.end_coord * (local_coord[1] + 1))
end

# Return a vector of coordinates the interpolation points of the geometry
function coords(a::Line2)
    return [a.start_coord, a.end_coord]
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::Line2)
    return 1
end

# Get the number of control points of an object
function number_of_control_points(a::Line2)
    return 2
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::Line2, position::Vector{Real})
    @assert(length(position) == 1, "Line2 is one dimensionsal.")
    return -1 <= position[1] <= 1
end

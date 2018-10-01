#===============================================================================
    PolyLine2.jl

    Represents a piecewise linear line in cartesian space

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

mutable struct PolyLine2 <: DiscreteGeometry3D
    coords :: Vector{Vector3D}

    function PolyLine2(a::Vector{Vector3D})
        new(a)
    end
end

# Return map a local coordinate to a point in space
function evaluate(a::PolyLine2, local_coord::Vector{Real})
    @assert(length(position) == 1, "PolyLine2 is one dimensionsal.")
    lower = floor(local_coord[1])
    upper = ceil(local_coord[1])
    return 0.5 * (a.coords[lower] * (local_coordlocal_coord[1] - 1) +
        a.coords[upper] * (local_coord[1] + 1))
end

# Return a vector of coordinates the interpolation points of the geometry
function coords(a::PolyLine2)
    return coords
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::PolyLine2)
    return 1
end

# Get the number of control points of an object
function number_of_control_points(a::PolyLine2)
    return length(coords)
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::PolyLine2, position::Vector{Real})
    @assert(length(position) == 1, "PolyLine2 is one dimensionsal.")
    return 0 <= position[1] <= length(a.coords)
end

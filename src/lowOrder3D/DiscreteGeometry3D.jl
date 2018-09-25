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

abstract type DiscreteGeometry3D
end

# Return map a local coordinate to a point in space
function evaluate(a::DiscreteGeometry3D, local_coord::Vector{Real})
    error("Not yet implemented for ", typeof(a), ".")
end

# Return a vector of coordinates the interpolation points of the geometry
function coords(a::DiscreteGeometry3D)
    return a.coords
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::DiscreteGeometry3D)
    error("Not yet implemented for ", typeof(a), ".")
end

# Get the number of control points of an object
function number_of_control_points(a::DiscreteGeometry3D)
    return length(a.coords)
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::DiscreteGeometry3D, position::Vector{Real})
    error("Not yet implemented for ", typeof(a), ".")
end

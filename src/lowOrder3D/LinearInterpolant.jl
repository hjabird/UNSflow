#===============================================================================
    LinearInterpolant.jl

    Represents a linear interpolation of two values in [-1,1]

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

mutable struct LinearInterpolant
    # Return map a local coordinate to a point in space
    function evaluate(
        a::LinearInterpolant,
        vals::Vector{T},
        local_coord::Vector{S}) where {T <: Any, S <: Real}

        @assert(length(local_coord) == 1, string("LinearInterpolant is 1"*
            " dimensionsal. Argument local_coord was ", length(local_coord),
            " long."))
        @assert(length(vals) == 1, string("LinearInterpolant is defined by"*
            " 2 points. Argument vals was ", length(vals), "long."))
        return 0.5 * (vals[1] * (local_coord[1] - 1) +
            vals[2] * (local_coord[1] + 1))
    end

    # Get the number of dimensions that the space operates in (ie, line->1, surf->2)
    function local_dimensions()
        return 1
    end

    # Get the number of control points of an object
    function number_of_control_points()
        return 2
    end

    # Test if a point is in the bounds defined by the object
    function in_bounds(position::Vector{T}) where T <: Real
        @assert(length(position) == 1, string("LinearInterpolant is 1"*
            " dimensionsal. Argument local_coord was ", length(local_coord),
            " long."))
        return -1 <= position[1] <= 1
    end
end

#===============================================================================
    EquationSurf.jl

    Represents a surface defined by expressions that take arguments in
    [-1, 1], [-1, 1]

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

mutable struct EquationSurf
    surf_eq :: Function

    function EquationSurf(
        x_def :: Function,
        y_def :: Function,
        z_def :: Function
        )

        msg = "EquationSurf definition functions must take Vector{Float64}"*
            " as an argument (Vector of length 2)."
        @assert(hasmethod(x_def, Tuple{Vector{Float64}}), msg)
        @assert(hasmethod(x_def, Tuple{Vector{Float64}}), msg)
        @assert(hasmethod(x_def, Tuple{Vector{Float64}}), msg)
        msg2 = "Equation surf definition functions must return something"*
            " that can be converted to a Float64."
        @assert(typeof(x_def([0.0,0.0])) >: Float64, msg2)
        @assert(typeof(y_def([0.0,0.0])) >: Float64, msg2)
        @assert(typeof(z_def([0.0,0.0])) >: Float64, msg2)
        function def3d(x::Vector{Float64})
            @assert(length(x) == 2, "EquationSurf is two dimensionsal.")
            ret = Vector3D(x_def(x), y_def(x), z_def(x))
            return ret
        end
        new(def3d)
    end

    function EquationSurf(
        def3d :: Function
        )
        msg = "EquationSurf definition functions must take Vector{Float64}"*
            " as an argument (Vector of length 2)."
        @assert(hasmethod(def3d, Tuple{Vector{Float64}}), msg)
        # It probably isn't too expensive to give their method a quick test,
        # right? Don't you miss just knonwing the return type.
        test_res = def3d([-0.5, 0.5])
        @assert(typeof(test_res) == Vector3D, "Expected the EquationSurf"*
            " equation definition to return a UNSflow.Vector3D")
        new(def3d)
    end

end

# Return map a local coordinate to a point in space
function evaluate(a::EquationSurf, local_coord::Vector{T}) where T <: Real
    @assert(length(local_coord) == 2, "EquationSurf is two dimensionsal.")
    return a.surf_eq(local_coord)
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::EquationSurf)
    return 2
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::EquationSurf, position::Vector{T}) where T <: Real
    @assert(length(position) == 2, "EquationSurf is two dimensionsal.")
    return (-1 <= position[1] <= 1) && (-1 <= position[2] <= 1)
end

function discretise(
    a::EquationSurf, ::Type{BilinearQuad},
    xgrid::Vector{Float64}, ygrid::Vector{Float64})

    @assert(all(-1 .<= xgrid .<= 1), "All grid coordinates must be in [-1,1]")
    @assert(all(-1 .<= ygrid .<= 1), "All grid coordinates must be in [-1,1]")
    @assert(length(xgrid) >= 2, "xgrid must have length >= 2")
    @assert(length(ygrid) >= 2, "ygrid must have length >= 2")
    @assert(issorted(xgrid) || issorted(xgrid, rev=true), "xgrid must be"*
        " sorted.")
    @assert(issorted(ygrid) || issorted(ygrid, rev=true), "ygrid must be"*
        " sorted.")

    n_output = Vector{BilinearQuad}(undef, (length(xgrid)-1) *
        (length(ygrid) - 1))
    for i = 1 : length(xgrid) - 1
        for j = 1 : length(ygrid) - 1
            idx = (i-1) * (length(ygrid)-1) + j
            c1 = evaluate(a, [xgrid[i], ygrid[j]])
            c2 = evaluate(a, [xgrid[i+1], ygrid[j]])
            c3 = evaluate(a, [xgrid[i+1], ygrid[j+1]])
            c4 = evaluate(a, [xgrid[i], ygrid[j+1]])
            n_output[idx] = BilinearQuad(c1, c2, c3, c4)
        end
    end
    for child in n_output
        tst = child.c1
        tst = child.c2
        tst = child.c3
        tst = child.c4
    end
    return n_output
end

function discretise(
    a::EquationSurf, ::Type{PolyLine2},
    path_waypoints::Matrix{Float64})

    @assert(size(path_waypoints)[2] == 2, "Expecting an num_waypoints by 2"*
        " matrix for path_waypoints.")
    @assert(size(path_waypoints)[1] >= 2, "There must be more than 1 waypoint.")
    points = Vector{Vector3D}(undef, size(path_waypoints)[1])
    for i = 1 : size(path_waypoints)[1]
        points[i] = evaluate(a, path_waypoints[i, :])
    end
    ret = PolyLine2(points)
    return ret
end

function discretise(
    a::EquationSurf, ::Type{Spline3D},
    path_waypoints::Matrix{Float64})

    @assert(size(path_waypoints)[2] == 2, "Expecting an num_waypoints by 2"*
        " matrix for path_waypoints.")
    @assert(size(path_waypoints)[1] >= 2, "There must be more than 1 waypoint.")
    points = Vector{Vector3D}(undef, size(path_waypoints)[1])
    for i = 1 : size(path_waypoints)[1]
        points[i] = evaluate(a, path_waypoints[i, :])
    end
    return Spline3D(points)
end

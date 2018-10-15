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
    # If not used, or not matching surf_eq these should be Nothing
    x_def :: Any
    y_def :: Any
    z_def :: Any

    function EquationSurf(
        x_def :: Function,
        y_def :: Function,
        z_def :: Function
        )

        msg = "EquationSurf definition functions must be able to take "*
            "Vector{Float64} as an argument (Vector of length 2). Idealy"*
            " it should take Vector{T} where T<:Real as an argument."
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
        new(def3d, x_def, y_def, z_def)
    end

    function EquationSurf(
        def3d :: Function
        )
        msg = "EquationSurf definition functions must take Vector{Float64}"*
            " as an argument (Vector of length 2). Ideally define per "*
            "compenent for Real."
        @assert(hasmethod(def3d, Tuple{Vector{Float64}}), msg)
        # It probably isn't too expensive to give their method a quick test,
        # right? Don't you miss just knonwing the return type.
        test_res = def3d([-0.5, 0.5])
        @assert(typeof(test_res) == Vector3D, string("Expected the "*
            "EquationSurf equation definition to return a UNSflow.Vector3D."*
            " Instead, a ", typeof(test_res), " was returned."))
        new(def3d, Nothing, Nothing, Nothing)
    end

end

# Return map a local coordinate to a point in space
function evaluate(
    a::EquationSurf,
    local_coord::Vector{T},
    boundscheck::Bool=true) where T <: Real

    @assert(length(local_coord) == 2, "EquationSurf is two dimensionsal.")
    if boundscheck
        @assert(in_bounds(a, local_coord), string("Trying to evaluate a point "*
            "outside the the bounds of a a EquationSurf object ([-1,1]"*
            ", [-1,1]). Value of argument for local coord given: ",
            local_coord))
    end
    ret_val = a.surf_eq(local_coord)
    if any(isnan.(ret_val)) @warn "Encountered NaN whilst evaluating"*
            " EquationSurf" end
    return ret_val
end

function evaluate(
    a::EquationSurf,
    local_coord::Tuple{T, T},
    boundscheck::Bool=true) where {T <: Real}
    evaluate(a, [local_coord[1], local_coord[2]], boundscheck)
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

function derivative(a::EquationSurf, direction::Int, local_coord::Vector{T}
    ; boundscheck::Bool=true, dt::Float64=0.0001) where T <: Real
    if boundscheck
        @assert(in_bounds(a, local_coord), string("Trying to evaluate a point "*
            "outside the the bounds of a a EquationSurf object ([-1,1]"*
            ", [-1,1]). Value of argument for local coord given: ",
            local_coord))
    end
    @assert(1 <= direction <= 2, string("The direction must be either 1 or 2,"*
        " referring to the indexes of the local_coord used to evaluate "*
        "EquationSurf. The given value was: ", direction))
    @assert(dt > 0, string("The finite difference spacing dt must be positive.",
        " Given: ", dt))
    @assert(dt < 1, string("The finite difference spacing dt must be small.",
        " Given value of ", dt, " is too big."))

    deriv = Vector3D(0,0,0)
    # First we can try automatic differentiation.
    if  all(map(x->hasmethod(x, Tuple{Vector{Real}}),
        (a.x_def, a.y_def, a.z_def)))
        if direction == 1
            deriv.x = ForwardDiff.derivative(x->a.x_def([x, local_coord[2]]),
                local_coord[1])
            deriv.y = ForwardDiff.derivative(x->a.y_def([x, local_coord[2]]),
                local_coord[1])
            deriv.z = ForwardDiff.derivative(x->a.z_def([x, local_coord[2]]),
                local_coord[1])
        else
            deriv.x = ForwardDiff.derivative(x->a.x_def([local_coord[1], x]),
                local_coord[2])
            deriv.y = ForwardDiff.derivative(x->a.y_def([local_coord[1], x]),
                local_coord[2])
            deriv.z = ForwardDiff.derivative(x->a.z_def([local_coord[1], x]),
                local_coord[2])
        end
    else
        # Dumb central differencing.
        dcoord = direction == 1 ? [dt, 0] : [0, dt]
        deriv = (evaluate(a, local_coord + dcoord) -
            evaluate(a, local_coord - dcoord)) / (2*dt)
    end
    return deriv
end

function normal(
    a::EquationSurf, local_coord::Vector{T};
    boundscheck::Bool=true) where T <: Real

    # We'll just get the tangents and then take the cross product.
    dir1 = derivative(a, 1, local_coord; boundscheck=boundscheck)
    dir2 = derivative(a, 2, local_coord; boundscheck=boundscheck)
    return unit(cross(dir1, dir2))
end

function discretise(
    a::EquationSurf, ::Type{BilinearQuad},
    xgrid::Vector{T}, ygrid::Vector{T}) where T <:Real
    
    r = discretise(a, BilinearQuadSurf, xgrid, ygrid)
    return convert(Matrix{BilinearQuad}, r)
end

function discretise(a::EquationSurf, ::Type{BilinearQuadSurf},
    xgrid::Vector{T}, ygrid::Vector{T}) where T <:Real
    @assert(all(-1 .<= xgrid .<= 1), "All grid coordinates must be in [-1,1]")
    @assert(all(-1 .<= ygrid .<= 1), "All grid coordinates must be in [-1,1]")
    @assert(length(xgrid) >= 2, "xgrid must have length >= 2")
    @assert(length(ygrid) >= 2, "ygrid must have length >= 2")
    @assert(issorted(xgrid) || issorted(xgrid, rev=true), "xgrid must be"*
        " sorted.")
    @assert(issorted(ygrid) || issorted(ygrid, rev=true), "ygrid must be"*
        " sorted.")

    coords = map(x->evaluate(a, x), [[x, y] for x in xgrid, y in ygrid])
    return BilinearQuadSurf(coords)
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

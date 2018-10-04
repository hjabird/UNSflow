#===============================================================================
    Spline3D.jl

    Represents a spline in 3D space. The splines in Dierckx put into 3D space.

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
import Dierckx
import QuadGK

mutable struct Spline3D
    x_spline :: Dierckx.Spline1D
    y_spline :: Dierckx.Spline1D
    z_spline :: Dierckx.Spline1D
    limits :: Tuple{Float64, Float64}
    definition_points :: Vector{Float64}

    function Spline3D(
        ts :: Vector{Float64}, xs :: Vector{Float64},
        ys :: Vector{Float64}, zs :: Vector{Float64};
        k::Int64=3, bc="extrapolate")

        @assert((length(ts) == length(xs)) && (length(ys) == length(zs)) &&
            (length(xs) == length(ys)),
            "All input vectors must be of the same length.")
        sx = Dierckx.Spline1D(ts, xs, k=3, bc="error")
        sy = Dierckx.Spline1D(ts, ys, k=3, bc="error")
        sz = Dierckx.Spline1D(ts, zs, k=3, bc="error")
        lims = extrema(ts)
        new(sx, sy, sz, lims, ts)
    end

    function Spline3D(
        xs :: Vector{Float64},
        ys :: Vector{Float64}, zs :: Vector{Float64};
        k::Int64=3, bc="extrapolate")

        @assert(length(xs) == length(ys) == length(zs),
            "All input vectors must be of the same length.")
        ts = collect(0.:length(xs) - 1)
        return Spline3D(ts, xs, ys, zs, k=k, bc=bc)
    end

    function Spline3D(
        points :: Vector{Vector3D};
        k::Int64=3, bc="extrapolate")
        xs = map(x->x.x, points)
        ys = map(x->x.y, points)
        zs = map(x->x.z, points)
        return Spline3D(xs, ys, zs, k=k, bc=bc)
    end
end

function (a::Spline3D)(t::Real)
    return Vector3D(a.x_spline(t), a.y_spline(t), a.z_spline(t))
end

function evaluate(a::Spline3D, t::Real)
    return a(t)
end

function (a::Spline3D)(t::Vector{T}) where T <: Real
    return map(x->evaluate(a,x), t)
end

function evaluate(a::Spline3D, t::Vector{T}) where T <: Real
    return a(t)
end

function derivative(a::Spline3D, t::Real, nu::Int64=1)
    return Vector3D(
        Dierckx.derivative(a.x_spline, t, nu=nu),
        Dierckx.derivative(a.y_spline, t, nu=nu),
        Dierckx.derivative(a.z_spline, t, nu=nu)
    )
end

function derivative(a::Spline3D, t::Vector{T}, nu::Int64=1) where T <: Real
    return map(x->derivative(a, x, nu=nu), t)
end

function spline_length(spl::Spline3D, a::Real, b::Real)
    @assert(a <= b, "Integration limits must be in order. Call with correct"*
        " order or use -integrate(a::Spline3D, b, a) if intentional.")
    @assert(spl.limits[1] <= a, "a must be morethan or equal to the lowest "*
        "point defined by the spline.")
    @assert(spl.limits[2] >= b, "b must be lessthan or equal to the greatest "*
        "point defined by the spline.")

    integrand = x->abs(derivative(spl, x))
    return QuadGK.quadgk(integrand, a, b)[1]
end

function distribute(
    spl::Spline3D,
    a::Real, b::Real,
    n::Int, endpoints::Bool=true, reltol=0.01)
    @assert(a <= b, "Integration limits must be in order. Call with correct"*
        " order or use -integrate(a::Spline3D, b, a) if intentional.")
    @assert(spl.limits[1] <= a, "a must be morethan or equal to the lowest "*
        "point defined by the spline.")
    @assert(spl.limits[2] >= b, "b must be lessthan or equal to the greatest "*
        "point defined by the spline.")
    @assert(n > 0, "There must be a positive number of points to distrubute.")
    @assert(!endpoints || n > 1, "If endpoints are included the number of"*
        " points to distribute must be more than 1.")

    len = spline_length(spl, a, b)
    spacing = len / (endpoints == true ? Float64(n-1) : Float64(n))
    points = Vector{Vector3D}(undef, n)
    guess = a
    target = (endpoints ? 0 : spacing * 0.5)
    for i = 1 : n
        iterc = 0
        old_guess = -999999
        while iterc < 30 && abs(old_guess - guess) / spacing > reltol
            old_guess = guess
            guess = guess - (spline_length(spl, a, guess) - target) /
                abs(derivative(spl, guess))
            guess = guess > spl.limits[2] ? spl.limits[2] : guess
            guess = guess < spl.limits[1] ? spl.limits[1] : guess
        end
        target += spacing
        points[i] = spl(guess)
    end
    return points
end

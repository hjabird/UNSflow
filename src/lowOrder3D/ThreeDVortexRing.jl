#===============================================================================
    ThreeDVortexRing.jl

    Representation of a singular vortex filament ring.

    Initial code: HJAB 2018

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
include("ThreeDVector.jl")
include("ThreeDVorticity.jl")

type ThreeDVortexRing <: ThreeDVorticity
    c1 :: ThreeDVector
    c2 :: ThreeDVector
    c3 :: ThreeDVector
    c4 :: ThreeDVector
    strength :: Float64

    function ThreeDVortexRing(
        corner1 :: ThreeDVector,
        corner2 :: ThreeDVector,
        corner3 :: ThreeDVector,
        corner4 :: ThreeDVector,
        strength :: Float64
        )
        return new(corner1, corner2, corner3, corne4, strength)
    end
end

function ThreeDVortexRing(
    corner1 :: ThreeDVector,
    corner2 :: ThreeDVector,
    corner3 :: ThreeDVector,
    corner4 :: ThreeDVector
    )
    return ThreeDVortexRing(corner1, corner2, corner3, corner4, 0.0)
end

function ThreeDVortexRing()
    c1 = ThreeDVector()
    c2 = ThreeDVector()
    c3 = ThreeDVector()
    c4 = ThreeDVector()
    return ThreeDVortexRing(c1, c2, c3, c4, 0.0)
end

function convert(
    ::Type{Vector{ThreeDStraightVortexFilament}},
    a::ThreeDVortexRing)
    b = Vector{ThreeDStraightVortexFilament}(4)
    b[1] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    b[2] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    b[3] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    b[4] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    return b
end

function centre(ring::ThreeDVortexRing)
    return (ring.c1 + ring.c2 + ring.c3 + ring.c4) / 4
end

function effective_radius(ring::ThreeDVortexRing)
    c = centre(ring)
    return maximum([
        abs(ring.c1 - c),
        abs(ring.c2 - c),
        abs(ring.c3 - c),
        abs(ring.c4 - c),
    ])
end

function vorticity(ring::ThreeDVortexRing)
    return mapreduce(vorticity, ThreeDVector(0,0,0), +,
        convert(Vector{ThreeDStraightVortexFilament}(ring)))
end

function induced_velocity(
    inducing_ring :: ThreeDVortexRing,
    measurement_loc :: ThreeDVector
    )
    fils = convert(Vector{ThreeDStraightVortexFilament}, ring)
    return mapreduce(x->induced_velocity(x, measurement_loc),
        ThreeDVector(0,0,0), +, fils)
end

function induced_velocity_curl(
    ring :: ThreeDVortexRing,
    measurement_point :: ThreeDVector
    )
    fils = convert(Vector{ThreeDStraightVortexFilament}, ring)
    return mapreduce(x->induced_velocity(x, measurement_loc),
        [0. 0. 0.; 0. 0. 0.; 0. 0. 0.], +, fils)
end

function euler!(a::ThreeDVortexRing, b::ThreeDVorticity, dt::Real)
    v1 = induced_velocity(b, a.c1)
    v2 = induced_velocity(b, a.c2)
    v3 = induced_velocity(b, a.c3)
    v4 = induced_velocity(b, a.c4)
    c1 += v1 * dt
    c2 += v2 * dt
    c3 += v3 * dt
    c4 += v4 * dt
    return
end

#= END ThreeDVortexRing ------------------------------------------------------=#

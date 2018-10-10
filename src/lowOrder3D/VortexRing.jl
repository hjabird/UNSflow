#===============================================================================
    VortexRing.jl

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

mutable struct VortexRing <: Vorticity3D
    geometry :: PolyLine2
    strength :: Float64

    function VortexRing(
        corner1 :: Vector3D,
        corner2 :: Vector3D,
        corner3 :: Vector3D,
        corner4 :: Vector3D,
        strength :: Float64
        )
        return new(PolyLine2([corner1, corner2, corner3, corner4, corner1]),
            strength)
    end

    function VortexRing(quad :: BilinearQuad, strength :: Float64)
        return new(PolyLine2([quad.c1, quad.c2, quad.c3, quad.c4, quad.c1]),
            strength)
    end
end

function VortexRing(quad :: BilinearQuad)
    return VortexRing(quad, 1.0)
end

function Base.convert(
    ::Type{Vector{StraightVortexFilament}},
    a::VortexRing)
    b = Vector{StraightVortexFilament}(4)
    coords = a.goemetry.coords
    b[1] = StraightVortexFilament(coords[1], coords[2], a.strength)
    b[2] = StraightVortexFilament(coords[2], coords[3], a.strength)
    b[3] = StraightVortexFilament(coords[3], coords[4], a.strength)
    b[4] = StraightVortexFilament(coords[4], coords[1], a.strength)
    return b
end

function centre(ring::VortexRing)
    coords = ring.goemetry.coords
    return (coords[1] + coords[2] + coords[3] + coords[4]) / 4
end

function normal(ring::VortexRing)
    g = ring.geometry.coords
    t1 = 0.25 * (g[4] + g[3] - g[2] - g[1])
    t2 = 0.25 * (g[2] + g[3] - g[1] - g[4])
    n = unit(cross(t1, t2))
    return n
end

function effective_radius(ring::VortexRing)
    c = centre(ring)
    coords = ring.goemetry.coords
    return maximum([
        abs(coords[1] - c),
        abs(coords[2] - c),
        abs(coords[3] - c),
        abs(coords[4] - c),
    ])
end

function vorticity(ring::VortexRing)
    return mapreduce(vorticity, Vector3D(0,0,0), +,
        convert(Vector{ThreeDStraightVortexFilament}(ring)))
end

function induced_velocity(
    inducing_ring :: VortexRing,
    measurement_loc :: Vector3D
    )
    fils = convert(Vector{ThreeDStraightVortexFilament}, ring)
    return mapreduce(x->induced_velocity(x, measurement_loc),
        Vector3D(0,0,0), +, fils)
end

function induced_velocity_curl(
    ring :: VortexRing,
    measurement_point :: Vector3D
    )
    fils = convert(Vector{StraightVortexFilament}, ring)
    return mapreduce(x->induced_velocity(x, measurement_loc),
        [0. 0. 0.; 0. 0. 0.; 0. 0. 0.], +, fils)
end

function euler!(a::VortexRing, b::Vorticity3D, dt::Real)
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

function state_vector_length(a::VortexRing)
    return 4 * 3
end

function vorticity_vector_length(this::VortexRing)
    return 1
end

function vorticity_vector(this::VortexRing)
    return [this.vorticity]
end

function update_using_vorticity_vector!(
    this::VortexRing,
    vort_vect::Vector{Float64})

    @assert(length(vort_vect) == 1)
    this.vorticity = vort_vect[1]
    return
end

function vorticity_vector_velocity_influence(
    this::VortexRing,
    mes_pnt::Vector3D
    )

    v = zeros(3, 1)
    old_vorticity = this.vorticity
    this.vorticity = 1.
    v = convert(Vector{Float}, induced_velocity(this, mes_pnt))
    this.vorticity = old_vorticity;
    return v
end
#= END VortexRing ------------------------------------------------------=#

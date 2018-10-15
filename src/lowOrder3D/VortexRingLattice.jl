#===============================================================================
    VortexRingLattice.jl

    A vortex ring lattice object - reduces computational overhead compared
    to discrete vortex rings. This is NOT a Vorticity3DCollector since it 
    doesn't have a children field.

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

mutable struct VortexRingLattice <: Vorticity3D
    strengths :: Matrix{Float64}
    geometry :: BilinearQuadSurf

    function VortexRingLattice(x_rings::Float64, y_rings::Float64)
        strs = zeros(x_rings, y_rings)
        geometry = BilinearQuadSurf(x_rings, y_rings)
        new(strs, geometry)
    end
    function VortexRingLattice(
        strengths :: Matrix{Float64},
        geometry :: BilinearQuadSurf)

        @assert(size(strengths) == size(geometry),
            string("Size of strengths != Size of geometry.",
            " size(strengths) = ", size(strengths),
            " and size(geometry) = ", size(geometry), "."))
        if any(isnan.(strengths)) || any(isinf.(strengths))
            @warn "strengths contained NaN or Inf value(s)."
        end
        new(strengths, geometry)
    end
end

#= VORTICITY3D functions -----------------------------------------------------=#
function centre(a::VortexRingLattice)
    c = coords(a.geometry)
    sum = mapreduce(x->x, +, c)
    sum /= length(c)
    return sum
end

function effective_radius(a::VortexRingLattice)
    middle = centre(a)
    rad = mapreduce(x->abs(middle - x), max, coords(a.geometry))
    return rad
end

function vorticity(a::VortexRingLattice)
    return Vector3D(0,0,0)
end

function induced_velocity(a::VortexRingLattice, measurement_loc::Vector3D)
    vel = Vector3D(0,0,0)
    id, jd = size(a.strengths)
    c = a.geometry.coordinates
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start = c[i, j]
            stop = c[i, j + 1]
            str = a.strengths[i - 1, j] - a.strengths[i, j]
            vel += induced_velocity(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.strengths[i, j] - a.strengths[i, j - 1]
            vel += induced_velocity(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # external horizontal filaments
    imax, jmax = size(a.ring_vertices)
    vel = mapreduce(
        i->induced_velocity(StraightVortexFilament, c[i, 1],
        c[i+1, 1], a.strengths[i, 1], measurement_loc),
        1 : id; init=vel)
    vel = mapreduce(
        i->induced_velocity(StraightVortexFilament, c[i, jmax],
        c[i+1, jmax], -a.strengths[i, jmax], measurement_loc),
        1 : id; init=vel)
    vel = mapreduce(
        j->induced_velocity(StraightVortexFilament, c[1, j],
        c[1, j+1], -a.strengths[1, j], measurement_loc),
        1 : jd; init=vel)
    vel = mapreduce(
        j->induced_velocity(StraightVortexFilament, c[imax, j],
        c[imax, j+1], a.strengths[imax, j], measurement_loc),
        1 : jd; init=vel)
    return vel
end

function induced_velocity_curl(a::VortexRingLattice)
    curl = zeros(3,3)
    id, jd = size(a.strengths)
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start =c[i, j]
            stop = c[i, j + 1]
            str = a.strengths[i - 1, j] - a.strengths[i, j]
            curl += induced_velocity_curl(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.strengths[i, j] - a.strengths[i, j - 1]
            curl += induced_velocity_curl(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # external horizontal filaments
    imax, jmax = size(c)
    curl = mapreduce(
        i->induced_velocity_curl(StraightVortexFilament, c[i, 1],
        c[i+1, 1], a.strengths[i, 1], measurement_loc),
        1 : id; init=curl)
    curl = mapreduce(
        i->induced_velocity_curl(StraightVortexFilament, c[i, jmax],
            c[i+1, jmax], -a.strengths[i, jmax], measurement_loc),
        1 : id; init=curl)
    curl = mapreduce(
        j->induced_velocity_curl(StraightVortexFilament, c[1, j],
        c[1, j+1], -a.strengths[1, j], measurement_loc),
        1 : jd; init=curl)
     curl = mapreduce(
        j->induced_velocity_curl(StraightVortexFilament, c[imax, j],
        c[imax, j+1], a.strengths[imax, j], measurement_loc),
         1 : jd; init=curl)
    return curl
end

function state_vector_length(a::VortexRingLattice)
    return 3 * length(a.geometry.coordinates)
end

function state_vector(a::VortexRingLattice)
    v = zeros(state_vector_length(a))
    for i = 1 : length(a.geometry.coordinates)
        v[(i-1) * 3 + 1 : i * 3] = convert(Vector{Float64}, 
            a.geometry.coordinates[i])
    end
    return v
end

function update_using_state_vector!(
    this::VortexRingLattice,
    state_vect::Vector{Float64})

    n = length(this.ring_vertices)
    @assert(length(state_vector) == n * 3, string(
        "State vector was incorrect length. Length was ", length(state_vect),
        " but should have been ", n * 3, "."))
    if any(isnan.(state_vect) || isinf.(state_vect))
        @warn "Infinite or NaN values in state vector."
    end
    for i = 1 :  n
        this.geometry.coordinates[i] = 
            convert(Vector3D, state_vect[i*3 - 2: i*3])
    end
    return
end

function state_time_derivative(
    this::VortexRingLattice,
    inducing_bodies::Vorticity3D)

    svtd = Vector{Float64}(undef, 3 * length(this.geometry.coordinates))
    for i = 1 : length(this.ring_vertices)
        svtd[i*3 - 2 : i*3] = convert(Vector{Float64},
             induced_velocity(this, inducing_bodies))
    end
    return svtd
end

function vorticity_vector_length(this::VortexRingLattice)
    return length(this.strengths)
end

function vorticity_vector(this::VortexRingLattice)
    return vec(this.strengths)
end

function update_using_vorticity_vector!(
    this::VortexRingLattice, vort_vect::Vector{Float64})

    @assert(length(vort_vect) == vorticity_vector_length(this),
        string("Incorrect vorticity_vector length for VortexRingLattice.",
        " Expected a length of ", vorticity_vector_length(this), " but ",
        " was given a vector of length ", length(vort_vect)))
    for i = 1 : length(vort_vect)
        this.strengths[i] = vort_vect[i]
    end
    return
end

function vorticity_vector_velocity_influence(
    this :: VortexRingLattice,
    mes_pnt :: Vector3D)

    error("Need to implement this...")
end

#= CONVERSION FUNCTIONS ------------------------------------------------------=#
function Base.convert(::Type{Matrix{VortexRing}}, a::VortexRingLattice)
    ret = Matrix{VortexRing}(undef, size(a.strengths))
    g = a.geometry.coordinates
    s = a.strengths
    ret = map(
        x->VortexRing(g[x[1]+1, x[2]], g[x[1]+1, x[2]+1], g[x[1], x[2]+1],
            g[x[1], x[2]], s[x[1], x[2]] ),
        [(i, j) for i = 1 : size(s,1), j = 1 : size(s, 2)]
        )
    return ret
end

function Base.convert(::Type{Vector{VortexRing}}, a::VortexRingLattice)
    return vec(Base.convert(Matrix{VortexRing}, a))
end

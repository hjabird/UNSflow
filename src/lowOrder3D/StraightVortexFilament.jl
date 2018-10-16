#===============================================================================
    StraightVortexFilament.jl

    Representation of a straight singular vortex filament.

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

mutable struct StraightVortexFilament <: Vorticity3D
    geometry :: Line2
    vorticity :: Float64

    function StraightVortexFilament(
        start_coord :: Vector3D,
        end_coord :: Vector3D,
        vorticity :: T
        ) where T <: Real
        new(Line2(start_coord, end_coord), Float64(vorticity))
    end
end

function StraightVortexFilament(
    start_coord :: Vector3D,
    end_coord :: Vector3D
    )
    return StraightVortexFilament(Line2(start_coord, end_coord), 0.0)
end

StraightVortexFilament() =
    StraightVortexFilament(
        Vector3D([0., 0., 0.]),
        Vector3D([1., 1., 1.]),
        0.0)

function to_particles(a::StraightVortexFilament, radius::Real,
    kernal::Vortex3DRegularisationFunctions)

    vec = a.geometry.end_coord - a.geometry.start_coord
    len = abs(vec)
    np = ceil(1.3 * len / (radius * 2))
    relative_div = len / np
    pos = map(
        x->a.geometry.start_coord + x * vec,
        relative_div .* collect(0.5 : 0.5 + np))
    vort = vec * a.vorticity / np
    ret = map(x->VortexParticle3D(x, vort, relative_div / 2, kernal), pos)
    return Vorticity3DSimpleCollector(ret)
end

function centre(filament::StraightVortexFilament)
    return (filament.geometry.start_coord + filament.geometry.end_coord) / 2
end

function effective_radius(filament::StraightVortexFilament)
    return (filament.geometry.start_coord - filament.geometry.end_coord) / 2
end

function vorticity(filament::StraightVortexFilament)
    return filament.vorticity * (filament.geometry.end_coord
        - filament.geometry.start_coord)
end

function induced_velocity(
    filament :: StraightVortexFilament,
    measurement_loc :: Vector3D
    )
    # Other objects - like vortex rings - will want to be able
    # to use this in functional way, hence the indirection.
    return induced_velocity(
        StraightVortexFilament, 
        filament.geometry.start_coord,
        filament.geometry.end_coord,
        filament.vorticity,
        measurement_loc
        )
end

function induced_velocity_curl(
    filament :: StraightVortexFilament,
    measurement_point :: Vector3D
    )
    # Other objects - like vortex rings - will want to be able
    # to use this in functional way, hence the indirection.
    return induced_velocity_curl(
        StraightVortexFilament, 
        filament.geometry.start_coord,
        filament.geometry.end_coord,
        filament.vorticity,
        measurement_loc
        )
end

function euler!(
    a::StraightVortexFilament,
    b::Vorticity3D,
    dt::Real)

    vels = induced_velocity(b, a.geometry.start_coord)
    vele = induced_velocity(b, a.geometry.end_coord)
    a.geometry.start_coord += vels * dt
    a.geometry.end_coord += vele * dt
    return
end

function state_vector_length(a::StraightVortexFilament)
    return 2 * 3
end

function vorticity_vector_length(this::StraightVortexFilament)
    return 1
end

function vorticity_vector(this::StraightVortexFilament)
    return [this.vorticity]
end

function update_using_vorticity_vector!(
    this::StraightVortexFilament,
    vort_vect::Vector{Float64})

    @assert(length(vort_vect) == 1)
    this.vorticity = vort_vect[1]
    return
end

function vorticity_vector_velocity_influence(
    this::StraightVortexFilament,
    mes_pnt::Vector3D
    )

    v = zeros(3, 1)
    old_vorticity = this.vorticity
    this.vorticity = 1.
    v = convert(Vector{Float}, induced_velocity(this, mes_pnt))
    this.vorticity = old_vorticity;
    return v
end

function induced_velocity(::Type{StraightVortexFilament},
        start::Vector3D, stop::Vector3D, strength::Float64,
        measurement_loc :: Vector3D)
    r0 = stop - start
    r1 = measurement_loc - start
    r2 = measurement_loc - stop
    if(abs(r1) <= eps(Float64) || abs(r2) <= eps(Float64))
        return Vector3D(0,0,0)
    end
    # From Katz & Plotkin, Eq(2.72), pg41
    term1 = strength / (4 * pi)
    term2n = cross(r1, r2)
    term2d = abs(cross(r1, r2)) ^ 2
    term3 = dot(r0, unit(r1) - unit(r2))
    vel =  term1 * (term2n / term2d) * term3
    vel = (!isfinite(vel) ? Vector3D(0,0,0) : vel)
    return vel
end

function induced_velocity_curl(
    ::Type{StraightVortexFilament},
    start::Vector3D, stop::Vector3D, strength::Float64,
    measurement_loc :: Vector3D
    )

    r0 = stop - start
    r1 = measurement_point - start
    r2 = measurement_point - stop
    # Notes, HJAB, Book 4, pg.42-pg.43 for derivation of the general theme
    # and pg70 for conversion to matrix expression.
    term1 = strength / ( 4 * pi)
    term211 = -r0 / (abs(cross(r1, r0))^2)
    term2121 = dot(r0, r1) / abs(r1)
    term2122 = -dot(r0, r2) / abs(r2)
    term221 = 3.0 / abs(r0)
    term2221 = abs(cross(r0, r1)) / abs(r1)
    term2222 = -abs(cross(r0, r2)) / abs(r2)
    A = term211 * (term2121 + term2122)
    B = term221 * (term2221 + term2222)
    C = term1 * [B -A.z A.y; A.z B -A.x; -A.y A.x B]
    return C
end
#= END StraightVortexFilament ------------------------------------------=#

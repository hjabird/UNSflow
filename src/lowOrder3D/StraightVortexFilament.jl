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
    return abs(filament.geometry.start_coord - filament.geometry.end_coord) / 2
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

function steady_force(a::StraightVortexFilament,
    vel_fn::Function, 
    density::Real=1, 
    samples::Int=1)

    return steady_force(StraightVortexFilament,
        a.geometry.start_coord,
        a.geometry.end_coord,
        a.vorticity,
        vel_fn,
        density,
        samples)
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
    #=if(abs(r1) <= eps(Float64) || abs(r2) <= eps(Float64))
        println("Here: ", r1, " ", r2)
        return Vector3D(0,0,0)
    end=#
    # From Katz & Plotkin, Eq(2.72), pg41
    term1 = strength / (4 * pi)
    term2n = cross(r1, r2)
    term2d = abs(cross(r1, r2)) ^ 2
    if term2d < 1e-8 
        # We might be on the filament or near to it - careful check time.
        d = abs(term2n) / abs(r0) 
        if d < eps(Float64) * abs(r0) * 10 || isnan(d)
            return Vector3D(0,0,0)
        end
    end
    term3 = dot(r0, unit(r1) - unit(r2))
    vel =  term1 * (term2n / term2d) * term3
    if !isfinite(vel)        
        d = abs(term2n) / abs(r0) 
        @warn string("Evaluated induced velocity from vortex filament as ",
            "non-finite. Filament vorticity density was ", vorticity, ", and",
            " distance of point from filament was ", d, ".")
    end
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

# The Kutta-Joukowski theorum
function steady_force(
    ::Type{StraightVortexFilament},
<<<<<<< HEAD
    start::Vector3D, stop::Vector3D, strength::Real, 
=======
    start::Vector3D, stop::Vector3D, strength::Float64, 
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
    vel_fn::Function, 
    density::Real=1, 
    samples::Int=1)

    @assert(hasmethod(vel_fn, Tuple{Vector3D}),
        "Expected induced_vel_fn passed to steady force function to take a "*
        "single arguement of type UNSflow.Vector3D.")
    @assert(samples > 0, string("Samples was expected to be 1 or more. Given ",
        "value was ", samples, "."))

    mps = collect(-1: 2/samples : 1)
    mes_locs = map(x->-0.5 * start * (x - 1) + 0.5 * stop * (x + 1), mps)
<<<<<<< HEAD
    mes_locs = (mes_locs[1:end-1] + mes_locs[2:end]) / 2
=======
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
    vels = map(vel_fn, mes_locs)
    force = mapreduce(x->density*cross(x, (stop - start) * strength / samples),
        +, vels, init=Vector3D(0,0,0))
    return force
end
#= END StraightVortexFilament ------------------------------------------=#

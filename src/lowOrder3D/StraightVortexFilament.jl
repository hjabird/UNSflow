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
    start_coord :: Vector3D
    end_coord :: Vector3D
    vorticity :: Float64

    function StraightVortexFilament(
        start_coord :: Vector3D,
        end_coord :: Vector3D,
        vorticity :: T
        ) where T <: Real
        new(start_coord, end_coord, Float64(vorticity))
    end
end

function StraightVortexFilament(
    start_coord :: Vector3D,
    end_coord :: Vector3D
    )
    return StraightVortexFilament(start_coord, end_coord, 0.0)
end

StraightVortexFilament() =
    StraightVortexFilament(
        Vector3D([0., 0., 0.]),
        Vector3D([1., 1., 1.]),
        0.0)

function centre(filament::StraightVortexFilament)
    return (filament.start_coord + filament.end_coord) / 2
end

function effective_radius(filament::StraightVortexFilament)
    return (filament.start_coord - filament.end_coord) / 2
end

function vorticity(filament::StraightVortexFilament)
    return filament.vorticity * (filament.end_coord - filament.start_coord)
end

function induced_velocity(
    filament :: StraightVortexFilament,
    measurement_loc :: Vector3D
    )
    r0 = filament.end_coord - filament.start_coord
    r1 = measurement_loc - filament.start_coord
    r2 = measurement_loc - filament.end_coord
    if(abs(r1) <= eps(Float64) || abs(r2) <= eps(Float64))
        return Vector3D(0,0,0)
    end
    # From Katz & Plotkin, Eq(2.72), pg41
    term1 = filament.vorticity / (4 * pi)
    term2n = cross(r1, r2)
    term2d = abs(cross(r1, r2)) ^ 2
    term3 = dot(r0, unit(r1) - unit(r2))
    vel =  term1 * (term2n / term2d) * term3
    vel = (!isfinite(vel) ? Vector3D(0,0,0) : vel)
    return vel
end

function induced_velocity_curl(
    filament :: StraightVortexFilament,
    measurement_point :: Vector3D
    )

    r0 = filament.end_coord - filament.start_coord
    r1 = measurement_point - filament.start_coord
    r2 = measurement_point - filament.end_coord
    # Notes, HJAB, Book 4, pg.42-pg.43 for derivation of the general theme
    # and pg70 for conversion to matrix expression.
    term1 = filament.vorticity / ( 4 * pi)
    term211 = -r0 / (abs(cross(r1, r0))^2)
    term2121 = dot(r0, r1) / abs(r1)
    term2122 = -dot(r0, r2) / abs(r2)
    term221 = 3.0 / abs(r0)
    term2221 = abs(cross(r0, r1)) / abs(r1)
    term2222 = -abs(cross(r0, r2)) / abs(r2)
    #term = term211 * (term2121 + term2122) +
    #    term221 * (term2221 + term2222)
    A = term211 * (term2121 + term2122)
    B = term221 * (term2221 + term2222)
    C = term1 * [B -A.z A.y; A.z B -A.x; -A.y A.x B]
    return C
end

function euler!(
    a::StraightVortexFilament,
    b::Vorticity3D,
    dt::Real)

    vels = induced_velocity(b, a.start_coord)
    vele = induced_velocity(b, a.end_coord)
    a.start_coord += vels * dt
    a.end_coord += vele * dt
    return
end

#= END StraightVortexFilament ------------------------------------------=#

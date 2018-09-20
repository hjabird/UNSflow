#===============================================================================
    ThreeDVortexParticle.jl

    Representation of a vortex particle.

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
include("ThreeDVector.jl")
include("ThreeDVortexRegularisationFunctions.jl")
include("ThreeDVorticity.jl")

type ThreeDVortexParticle <: ThreeDVorticity
    coord :: ThreeDVector
    vorticity :: ThreeDVector
    radius :: Float64

    kernel_funcs :: ThreeDVortexRegularisationFunctions

    function ThreeDVortexParticle(
        coordinate::ThreeDVector,
        vorticity_vector::ThreeDVector,
        particle_radius::Real,
        kernel_functions=threed_winckelmans_kernels()
        )
        @assert(particle_radius >= 0.0)
        new(coordinate, vorticity_vector, particle_radius, kernel_functions)
    end
end

function centre(a::ThreeDVortexParticle)
    return coord
end

function effective_radius(a::ThreeDVortexParticle)
    return a.radius * a.kernel_funcs.radius_modifier()
end

function vorticity(a::ThreeDVortexParticle)
    return a.vorticity
end

function induced_velocity(
    particle::ThreeDVortexParticle,
    measurement_point::ThreeDVector)
    # Acceleration method & avoid singularities which most people are probably
    # evaluating accidently anyway.
    if iszero(particle.vorticity) || particle.coord == measurement_point
        return ThreeDVector(0, 0, 0)
    end
    rad = particle.coord - measurement_point
    a = particle.kernel_funcs.g(abs(rad)/particle.radius) / (4. * pi)
    den = abs(rad)^3
    c = cross(rad, particle.vorticity)
    vel = c * a / den
    return vel
end

function induced_velocity_curl(
    particle::ThreeDVortexParticle,
    measurement_point::ThreeDVector)

    # This is taken from Robertson and Joo, 2010 "Vortex Particle Aerodynamic
    # Modelling of Perching Manoeuvres with Micro Air Vehicles", and I've
    # factorised out the term for the affected particles vorticity, such that
    # it turns into a matrix multiplication dalpha/dt = A x alpha where
    # A is what we are generating here.
    if iszero(particle.vorticity) || particle.coord == measurement_point
        return zeros(3,3)
    end

    rad = particle.coord - measurement_point
    rho = abs(rad)/particle.radius
    gf = particle.kernel_funcs.g(rho)
    zetaf = particle.kernel_funcs.zeta(rho)
    omega = particle.vorticity
    sigma = particle.radius
    r_cross_k = cross(rad, omega)

    t1 = 1. / (4 * pi * sigma ^ 3)
    t21n = omega * gf   # Take out the negative since cross is anticommutative
    t21d = (abs(rad) / sigma) ^ 3
    t21 = t21n / t21d
    t21mat = [0.0 -t21.z t21.y; t21.z 0.0 -t21.x; -t21.y t21.x 0.0]

    t221 = 1 / (abs(rad)^2)
    t222 = 3 * gf / (abs(rad) / sigma) ^ 3 - zetaf
    t223 = outer(r_cross_k, rad)
    t22 = t221 * t222 * t223

    # Dot product: a.b -> [a.x,0,0; 0,a.y,0; 0,0,a.z]{b}
    # cross product: axb -> [0, -a.z, a.y; a.z,0,-a.x; -a.y,a.x,0 ]{b}
    A = t1 * (t22 + t21mat)
    return A
end

function euler!(a::ThreeDVortexParticle, b::ThreeDVorticity, dt::Real)
    vel = induced_velocity(b, a.coord)
    dvort = induced_velocity_curl(b, a.coord) * a.vorticity
    a.coord += vel * dt
    a.vorticity += dvort * dt
    return
end

#= END ThreeDVortexParticle --------------------------------------------------=#

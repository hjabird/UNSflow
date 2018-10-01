#===============================================================================
    VortexParticle3D.jl

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

mutable struct VortexParticle3D <: Vorticity3D
    geometry :: Point3D
    vorticity :: Vector3D
    radius :: Float64

    kernel_funcs :: Vortex3DRegularisationFunctions

    function VortexParticle3D(
        coordinate::Vector3D,
        vorticity_vector::Vector3D,
        particle_radius::Real,
        kernel_functions=threed_winckelmans_kernels()
        )
        @assert(particle_radius >= 0.0)
        new(Point3D(coordinate), vorticity_vector,
            particle_radius, kernel_functions)
    end
end

function centre(a::VortexParticle3D)
    return geometry.coord
end

function effective_radius(a::VortexParticle3D)
    return a.radius * a.kernel_funcs.radius_modifier()
end

function vorticity(a::VortexParticle3D)
    return a.vorticity
end

function induced_velocity(
    particle::VortexParticle3D,
    measurement_point::Vector3D)
    # Acceleration method & avoid singularities which most people are probably
    # evaluating accidently anyway.
    if iszero(particle.vorticity) || particle.geometry.coord==measurement_point
        return Vector3D(0, 0, 0)
    end
    rad = particle.geometry.coord - measurement_point
    a = particle.kernel_funcs.g(abs(rad)/particle.radius) / (4. * pi)
    den = abs(rad)^3
    c = cross(rad, particle.vorticity)
    vel = c * a / den
    return vel
end

function induced_velocity_curl(
    particle::VortexParticle3D,
    measurement_point::Vector3D)

    # This is taken from Robertson and Joo, 2010 "Vortex Particle Aerodynamic
    # Modelling of Perching Manoeuvres with Micro Air Vehicles", and I've
    # factorised out the term for the affected particles vorticity, such that
    # it turns into a matrix multiplication dalpha/dt = A x alpha where
    # A is what we are generating here.
    if iszero(particle.vorticity) || particle.geometry.coord==measurement_point
        return zeros(3,3)
    end

    rad = particle.geometry.coord - measurement_point
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

function euler!(a::VortexParticle3D, b::Vorticity3D, dt::Real)
    vel = induced_velocity(b, a.geometry.coord)
    dvort = induced_velocity_curl(b, a.geometry.coord) * a.vorticity
    a.geometry.coord += vel * dt
    a.vorticity += dvort * dt
    return
end

function state_vector(this::VortexParticle3D)
    # Because of the way vortex particles work, we it has to be redefined
    # there too
    state_vect = Vector{Float64}(undef, 6)
    state_vect[1:3] = convert(Vector{Float64}, this.geometry.coord)
    state_vect[4:6] = convert(Vector{Float64}, this.vorticity)
    return state_vect
end

function update_using_state_vector(
    this::VortexParticle3D,
    state_vect::Vector{Float64})

    state_vect[1:3] = convert(Vector{Float64}, this.geometry.coord)
    state_vect[4:6] = convert(Vector{Float64}, this.vorticity)
    state_vect = state_vect[7:end]
    return state_vect
end

function state_time_derivative(
    this::VortexParticle3D,
    inducing_bodies::Vorticity3D)
    # We assume here that the change in due only due to convection of points.
    # Otherwise we need to specially define this such as for VortexParticle3D
    deriv_vect = Vector{Float64}(undef, 6)
    deriv_vect[1:3] = convert(Vector{Float64},
            induced_velocity(inducing_bodies, this.geometry.coord))
    deriv_vect[4:6] = convert(Vector{Float64},
            induced_velocity_curl(inducing_bodies, this.geometry.coord) *
            this.vorticity )
    return deriv_vect
end

#= END VortexParticle3D --------------------------------------------------=#

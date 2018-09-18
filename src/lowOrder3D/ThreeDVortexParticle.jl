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

type ThreeDVortexParticle <: ThreeDVorticityBody
    coord :: ThreeDVector
    vorticity :: ThreeDVector
    radius :: Float64

    kernel_funcs :: ThreeDVortexParticleKernelFunctions

    velocity :: ThreeDVector
    vorticity_time_derivative :: ThreeDVector
end

function centre(a::ThreeDVortexParticle)
    return coord
end

function effective_radius(a::ThreeDVortexParticle)
    return a.radius * a.kernel_funcs.radius_modifier()
end

function induced_velocity(particle::ThreeDVortexParticle,
    measurement_point::ThreeDVector)
    # Acceleration method & avoid singularities which most people are probably
    # evaluating accidently anyway.
    if iszero(particle.vorticity) || a.coord == measurement_point
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

    rad = particle.coord - measurement_point
    rho = abs(rad)/particle.radius
    gf = particle.kernel_funcs.g(rho)
    zetaf = particle.kernel_funcs.zeta(rho)
    omega = particle.vorticity
    r_cross_k = cross(rad, omega)

    t1 = 1. / (4 * pi * sigma_k ^ 3)
    t21n = g * vorticity  # Take out the negative since cross is anticommutative
    t21d = (abs(r) / sigma_k) ^ 3
    t21 = t21n / t21d

    t221 = 1 / (abs(r)^2)
    t222 = 3 * g / (abs(r) / sigma_k) ^ 3 - f
    t223 = r * r_cross_k
    t22 = t221 * t222 * t223
    # Dot product: a.b -> [a.x,0,0; 0,a.y,0; 0,0,a.z]{b}
    # cross product: axb -> [0, -a.z, a.y; a.z,0,-a.x; -a.y,a.x,0 ]{b}
    A = t1 * [t22.x -t21.z t21.y; t21.z t22.y -t21.x; -t21.y t21.x t22.z]
    return A
end


"""
According to input vortex particle field particles, and kernal functions g and
f, compute the change for a timestep dt and return this as a new
particles_updated array of particles.
"""
function one_step(
    particles::Array{ThreeDVortexParticle},
    dt::Float64,
    g_function::Function,
    f_function::Function
    )

    delta_x = Array{ThreeDVector, 1}(size(particles))
    delta_vort = Array{ThreeDVector, 1}(size(particles))
    particles_updated = deepcopy(particles)

    function get_dx(
        particle_idx::Int64,
        particles::Array{ThreeDVortexParticle})

        location = particles[particle_idx].coord
        v = map(x->induced_velocity(x, location, g_function),
                        particles[vcat(1:particle_idx-1, particle_idx+1:end)])
        dx = sum(v) * dt
        return dx
    end

    function get_dvort(
        particle_idx::Int64,
        particles::Array{ThreeDVortexParticle})

        location = particles[particle_idx].coord
        vo = map(x->rate_of_change_of_vorticity(
                    x, particles[particle_idx], g_function, f_function),
                    particles[vcat(1:particle_idx-1, particle_idx+1:end)])
        dvo = sum(vo) * dt
        return dvo
    end

    delta_x = map(i->get_dx(i, particles), 1:length(particles))
    delta_vort = map(i->get_dvort(i, particles), 1:length(particles))

    for i = 1 : length(particles)
        particles_updated[i].coord += delta_x[i]
        particles_updated[i].vorticity += delta_vort[i]
    end
    return particles_updated
end

"""
Compute the velocity induced by ThreeDParticle particle at ThreeDCoordinate
coordinate, given a g function.
"""
function induced_velocity(
    particle::ThreeDVortexParticle,
    coordinate::ThreeDVector,
    g_function::Function
    )
    # Robertson 2010 Eq. 3
    if iszero(particle.vorticity)
        return ThreeDVector(0, 0, 0)
    end
    rad = particle.coord - coordinate
    a = g_function(abs(rad)/particle.radius) / (4. * pi)
    den = abs(rad)^3
    c = cross(rad, particle.vorticity)
    vel = c * a / den
    return vel
end

"""
The vorticity of a particle over time changes as vortices stretch
 and what not. This RETURNS (doesn't change the value of) the rate of change
 of particle j with respect to time as domega_x / dt, domega_y / dt,
 domega_z / dt
"""
function rate_of_change_of_vorticity(
    particle_j::ThreeDVortexParticle,
    particle_k::ThreeDVortexParticle,
    g_function::Function,
    f_function::Function
    )
    if iszero(particle_j.vorticity) || iszero(particle_k.vorticity)
        return ThreeDVector(0, 0, 0)
    end
    # Robertson 2010 Eq. 4
    sigma_k = particle_k.radius
    om_j = particle_j.vorticity
    om_k = particle_k.vorticity
    r = particle_j.coord - particle_k.coord
    g = g_function(abs(r) / particle_k.radius)
    f = f_function(abs(r) / particle_k.radius)
    om_cross = cross(particle_j.vorticity, particle_k.vorticity)
    r_cross_k = cross(r, particle_k.vorticity)
    om_dot_r = dot(particle_j.vorticity, r)

    t1 = 1. / (4 * pi * sigma_k ^ 3)
    t21n = - g * om_cross
    t21d = (abs(r) / sigma_k) ^ 3
    t21 = t21n / t21d
    t221 = 1 / (abs(r)^2)
    t222 = 3 * g / (abs(r) / sigma_k) ^ 3 - f
    t223 = om_dot_r * r_cross_k
    t22 = t221 * t222 * t223
    t2 = t21 + t22
    t = t1 * t2
    return t
end

"""Take a set of vortex particles and the value of g for each of them,
 and computes the velocity of each. Returns a vector of velocities. """
function particle_velocities(
    particles::Array{ThreeDVortexParticle},
    particle_g_function::Function)

    @assert(particles.size() == particles_g_values.size())  # Same size arrays

    vel = zeros(particles.size())
    i :: Int64
    v0 :: ThreeDVector
    for i = 1 : particles.size()
        x0 = particles[i].coord
        zero!(v0)
        vel[i] = mapreduce(x->induced_velocity(x, x0, particle_g_function),
            +, v0, particles)
    end
    return vel
end

"""Take a set of vortex particles and the value of g and f for each of them,
 and computes the rate of change of vorticity for each of them. """
function particle_rate_of_change_of_vorticity(
    particles::Array{ThreeDVortexParticle},
    particle_g_function::Function,
    particle_f_function::Function)

    @assert(particles.size() == particles_g_values.size())  # Same size arrays

    dvort = zeros(particles.size())
    i :: Int64
    v0 :: ThreeDVector
    for i = 1 : particles.size()
        x0 = particles[i].coord
        zero!(v0)
        dvort[i] = mapreduce(
            x->rate_of_change_of_vorticity(dvort[i], x[1],
                                        particle_g_values, particle_f_values),
            +, v0, particles)
    end
    return vel
end
#= END ThreeDVortexParticle --------------------------------------------------=#

#===============================================================================
    ThreeDVortexParticle methods
------------------------------------------------------------------------------=#
"""The induced velocity at a point 'coordinate' due to a vortex particle.
reduction_factor_fn is the function g that avoids velocity singularities."""
function ind_vel(
    particle::ThreeDVortexParticle,
    coordinate::ThreeDVector,
    reduction_factor_fn::Function
     )
    if particle.coord == coordinate || abs(particle.vorticity) == 0.0
    	return ThreeDVector(0., 0., 0.)
    end
    # Robertson 2010 Eq. 3 with additional - sign.
    rad = particle.coord - coordinate
    a = reduction_factor_fn(abs(rad)/particle.size) / (4. * pi)
    den = abs(rad)^3
    c = cross(rad, particle.vorticity)
    vel :: ThreeDVector = c * a / den
    return vel
end

function ind_vel(
    particles :: Vector{ThreeDVortexParticle},
    coordinate :: ThreeDVector,
    reduction_factor_fn :: Function
     )
     z = ThreeDVector(0., 0., 0.)
     v = mapreduce(x->ind_vel(x, coordinate, reduction_factor_fn),
        +, z, particles)
     return v
 end

 function ind_vel(
     particles :: Vector{ThreeDVortexParticle},
     reduction_factor_fn :: Function
      )
      v = map(x->ind_vel(particles, x.coord, reduction_factor_fn), particles)
      return v
  end

""" The time rate of change of vorticity of a particle
"dvort_induced_on_particle" due to another particle "inducing particle" """
function ind_dvortdt(
    dvort_induced_on_particle :: ThreeDVortexParticle,
    inducing_particle :: ThreeDVortexParticle,
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )
    if dvort_induced_on_particle.coord == inducing_particle.coord ||
	abs(inducing_particle.vorticity) == 0.0
        return ThreeDVector(0., 0., 0.)
    end

    # Robertson 2010 Eq. 4
    sigma_k = inducing_particle.size
    om_j = dvort_induced_on_particle.vorticity
    om_k = inducing_particle.vorticity
    r = dvort_induced_on_particle.coord - inducing_particle.coord
    g = reduction_factor_fn(abs(r) / inducing_particle.size)
    f = vorticity_fraction_fn(abs(r) / inducing_particle.size)
    om_cross = cross(dvort_induced_on_particle.vorticity,
                            inducing_particle.vorticity)
    r_cross_k = cross(r, inducing_particle.vorticity)
    om_dot_r = dot(dvort_induced_on_particle.vorticity, r)

    t1 = 1. / (4 * pi * sigma_k ^ 3)
    t21n = - g * om_cross
    t21d = (abs(r) / sigma_k) ^ 3
    t21 = t21n / t21d
    t221 = 1 / (abs(r)^2)
    t222 = 3 * g / (abs(r) / sigma_k) ^ 3 - f
    t223 = om_dot_r * r_cross_k
    t22 = t221 * t222 * t223
    t2 = t21 + t22
    t = t1 * t2
    return t
end

function ind_dvortdt(
    dvort_induced_on_particle :: ThreeDVortexParticle,
    inducing_particles :: Vector{ThreeDVortexParticle},
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )
    z = ThreeDVector(0., 0., 0.)
    v = mapreduce(x->ind_dvortdt(dvort_induced_on_particle, x,
        reduction_factor_fn, vorticity_fraction_fn), +, z, inducing_particles)
    return v
end

function ind_dvortdt(
    particles :: Vector{ThreeDVortexParticle},
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )
    v = map(x->ind_dvortdt(x, particles,
        reduction_factor_fn, vorticity_fraction_fn), particles)
    return v
end

"""
Compute the velocity and rate of change of vorticity of a set of vortex
particles on itself.
"""
function mutual_ind(
    particles::Vector{ThreeDVortexParticle},
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )

    vel = ind_vel(particles, reduction_factor_fn)
    dvortdt = ind_dvortdt(particles, reduction_factor_fn, vorticity_fraction_fn)
    for (particle, v, dvort) in zip(particles, vel, dvortdt)
        particle.velocity = v
        particle.vorticity_time_derivative = dvort
    end
    return particles
end

"""
Integrate the particle velocities and rate of change of vorticities over time
dt to compute new coordinate and vorticity. Assumes derivatives are already
computed using mutual_ind.
"""
function integration_step!(
    particles :: Vector{ThreeDVortexParticle},
    delta_t :: Float64
    )
    for particle in particles
        particle.coord += particle.velocity * dt
        particle.vorticity += particle.vorticity_time_derivative * dt
    end
    return particles
end

"""
Use the forward Euler method to calculate the location of vortex particles
at time dt in the particles future. Updates particles.
"""
function euler_forward_step!(
    particles :: Vector{ThreeDVortexParticle},
    delta_t :: Float64,
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )

    mutual_ind(particles, reduction_factor_fn, vorticity_fraction_fn)
    integration_step!(particles, dt)
    return particles
end

"""
Use the explicit midpoint method to calculate the location of vortex particles
at time dt in the particles future. Updates particles.
"""
function explicit_midpoint_step!(
    particles :: Vector{ThreeDVortexParticle},
    delta_t :: Float64,
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )

    mutual_ind(particles, reduction_factor_fn, vorticity_fraction_fn)
    particles_cpy = deepcopy(particles)
    integration_step!(particles_cpy, dt / 2.)
    mutual_ind(particles_cpy, reduction_factor_fn, vorticity_fraction_fn)
    for (particle, particlecpy) in zip(particles, particles_cpy)
        particle.velocity = particlecpy.velocity
        particle.vorticity_time_derivative =
            particlecpy.vorticity_time_derivative
    end
    integration_step!(particles, dt)
    return particles
end
 # END ThreeDVortexParticle methods ===========================================#

#===============================================================================
VortexParticle3DFlowFeatures.jl

Generate initial conditions for vortex particle test cases.

HJA Bird 2018
h.bird.1@research.gla.ac.uk

- vortex_particle_ring
- vortex_particle_sheet

===============================================================================#

""" Generates an array of vortex particles representing a vortex ring"""
function vortex_particle_ring(
    centre::Vector3D,
    normal_to_ring::Vector3D,
    radius::Float64,
    strength::Float64,
    num_particles::Int64,
    kernel::Vortex3DRegularisationFunctions
    )

    @assert(radius > 0)
    @assert(abs(normal_to_ring) > 0)
    @assert(num_particles >= 3)
    particles = Array{VortexParticle3D, 1}(undef, num_particles)

    # Angular definition
    dtheta = 2 * pi / num_particles
    theta_pos = dtheta * (1 : num_particles)
    particle_strength = radius * strength * dtheta
     # Close enough to the their segment of the circumference.
    particle_size = radius * dtheta * 1.3
    # Generate in plane coordinates
    coords = zeros(2, num_particles)
    for i = 1 : num_particles
        coords[:, i] = radius * [cos(theta_pos[i]), sin(theta_pos[i])]
    end

    # Generate a plane normal to n described by two orthogonal unit vectors
    # a and b.
    n = unit(normal_to_ring)
    if n[3] != 0
        a = Vector3D(1, 1, (n[1] + n[2])/n[3])
    elseif n[2] != 0
        a = Vector3D(1, (n[1] + n[3])/n[2], 1)
    elseif n[1] != 0
        a = Vector3D((n[2] + n[3])/n[1], 1, 1)
    end
    unit!(a)
    b = unit(cross(a, n))
    # Place particles in the new plane
    for i = 1 : num_particles
        from_centre = Vector3D(coords[1, i] * a + coords[2, i] * b)
        p = VortexParticle3D(
            centre + from_centre, # coord
            particle_strength * unit(cross(from_centre, n)), # vorticity
            particle_size, # size
            kernel
            )
        particles[i] = p
    end
    particle_collector = Vorticity3DSimpleCollector(particles)
    return particle_collector
end

"""Generates a sheet of vortex particles

Technically wrong with vortex particle strenght for now.

A function geometry_func f(x,y) returns Vector3D locations for the sheet
within the bounds given. The bounds array is intreted as of the
for [minx, maxx, miny, maxy]. Similalry, the vorticity density describes is of
the form f(x, y) and returns a Vector3D. num_x and num_y describe the number
of vortex particles to use. Returns Array{VortexParticle3D, 1}
"""
function vortex_particle_sheet(
    geometry_func::Function, # func(x,y) returns Vector3D
    vorticity_func::Function, # func(x,y) returns Vector3D
    bounds::Array{Float64}, # [minx, maxx, miny, maxy]
    num_x::Int64,
    num_y::Int64,
    kernel::Vortex3DRegularisationFunctions
    )

    @assert(num_x > 0)
    @assert(num_y > 0)
    @assert(bounds[1] < bounds[2])
    @assert(bounds[3] < bounds[4])
    linspace = (b1, b2, n) -> Vector{Float64}(b1 : (b2-b1)/n : b2)
    xs = linspace(bounds[1], bounds[2], num_x)
    ys = linspace(bounds[3], bounds[4], num_y)
    num_particles = num_x * num_y
    particles = Array{VortexParticle3D, 1}(undef, num_particles)

    for i = 1:num_x
        for j = 1:num_y
            z = Vector3D(0., 0., 0.)
            coord = geometry_func(xs[i], ys[j])
            # So this will be wrong for now.
            patch_size = ((bounds[2]-bounds[1])/num_x) *
                ((bounds[4]-bounds[3])/num_y)
            vort = vorticity_func(xs[i], ys[j]) * patch_size
            vort_size = 1.3 * sqrt(patch_size) / 2
            particles[num_x * (j - 1) + i] = VortexParticle3D(
                coord, vort, vort_size, kernel)
        end
    end

    return Vorticity3DSimpleCollector(particles)
end

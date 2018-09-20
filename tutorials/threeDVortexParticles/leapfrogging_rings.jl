#===============================================================================
leapfrogging_rings.jl

A demo of leapfrogging vortex rings using vortex particles.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
include("../../src/lowOrder3D/ThreeDVorticitySimpleCollector.jl")
include("../../src/lowOrder3D/ThreeDVortexParticle.jl")
include("ThreeDVortexParticleFlowFeatures.jl")
using WriteVTK  # We'll use this package to output to VTK for visualisation.

#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
num_steps = 2000
dt = 0.01
# Data saving parameters
basepath = "./output/leapfrogging_rings_"   # Where to write our output files
save_every = 10                       # Save every 10 steps.
# Initial conditions
ring_strength = [1., 1.]
ring_particles = [30, 30]
ring_radii = [1., 1.]
ring_locations = [0., 1.]

#=---------------------- Automated problem setup -----------------------------=#
num_rings = size(ring_particles)[1]
particles = ThreeDVorticitySimpleCollector()
for i = 1 : num_rings
    c = ThreeDVector(ring_locations[i], 0., 0.)
    n = ThreeDVector(1., 0., 0.)
    r = ring_radii[i]
    strength = ring_strength[i]
    n_ring_particles = ring_particles[i]
    ring = vortex_particle_ring(
        c, n, r, strength, n_ring_particles, threed_winckelmans_kernels()
        )
    particles = ThreeDVorticitySimpleCollector(
        get_children(particles),
        get_children(ring))
end
num_particles = size(particles)[1]

#=-------------------------- ODE time integration ----------------------------=#
println("Integrating in time...")
tic()
for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        points = zeros(3, num_particles)
        point_vorticity = zeros(3, num_particles)
        cells = Array{MeshCell, 1}(num_particles)
        for j = 1 : num_particles
            points[:,j] = [particles[j].coord.x,
                particles[j].coord.y, particles[j].coord.z]
            point_vorticity[:, j] = [particles[j].vorticity.x,
                particles[j].vorticity.y, particles[j].vorticity.z]
            cells[j] = MeshCell(VTKCellTypes.VTK_VERTEX, [j])
        end
        vtkfile = vtk_grid(string(basepath, i), points, cells)
        vtk_point_data(vtkfile, point_vorticity, "vorticity")
        outfiles = vtk_save(vtkfile)
    end

    # Calculate the next iteration
    euler!(particles, particles, dt)
end
toc()
#=------------------- Now repeat until it doesn't blow up --------------------=#

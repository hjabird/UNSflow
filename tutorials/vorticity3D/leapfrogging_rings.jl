#===============================================================================
leapfrogging_rings.jl

A demo of leapfrogging vortex rings using vortex particles.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
push!(LOAD_PATH, "../../src/")
import UNSflow
include("VortexFlowFeatures.jl")

let # The change in scope rules in Julia_1.0 make for a pickle otherwise.
#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
num_steps = 100
dt = 0.01
# Data saving parameters
basepath = "./output/leapfrogging_rings_"   # Where to write our output files
save_every = 10                       # Save every 10 steps.
# Initial conditions
ring_strength = [1., 1.]
ring_particles = [50, 50]
ring_radii = [1., 1.]
ring_locations = [0., 1.]

#=---------------------- Automated problem setup -----------------------------=#
num_rings = length(ring_particles)
particles = UNSflow.Vorticity3DSimpleCollector()
for i = 1 : num_rings
    c = UNSflow.Vector3D(ring_locations[i], 0., 0.)
    n = UNSflow.Vector3D(1., 0., 0.)
    r = ring_radii[i]
    strength = ring_strength[i]
    n_ring_particles = ring_particles[i]
    ring = vortex_particle_ring(
        c, n, r, strength, n_ring_particles, 
        UNSflow.threed_winckelmans_kernels())
    particles = UNSflow.Vorticity3DSimpleCollector(
        UNSflow.get_children(particles),
        UNSflow.get_children(ring))
end
num_particles = length(particles)

#=-------------------------- ODE time integration ----------------------------=#
@time for i = 1 : num_steps
    # Save only on some steps.
    if (i - 1) % save_every == 0
        mesh = UNSflow.UnstructuredMesh()
        push!(mesh, particles)
        UNSflow.to_vtk_file(mesh, string("output/leapfrogging_rings_", i))
    end

    # Calculate the next iteration
    UNSflow.euler!(particles, particles, dt)
end
#=------------------- Now repeat until it doesn't blow up --------------------=#
end

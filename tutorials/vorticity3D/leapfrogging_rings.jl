#===============================================================================
leapfrogging_rings.jl

A demo of leapfrogging vortex rings using vortex particles.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
include("../../src/lowOrder3D/Vorticity3DSimpleCollector.jl")
include("../../src/lowOrder3D/VortexParticle3D.jl")
include("../../src/lowOrder3D/DiscreteGeometry3DToVTK.jl")
include("VortexFlowFeatures.jl")
import WriteVTK  # We'll use this package to output to VTK for visualisation.

let # The change in scope rules in Julia_1.0 make for a pickle otherwise.
#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
num_steps = 2000
dt = 0.01
# Data saving parameters
basepath = "./output/leapfrogging_rings_"   # Where to write our output files
save_every = 10                       # Save every 10 steps.
# Initial conditions
ring_strength = [1., 1.]
ring_particles = [5, 5]
ring_radii = [1., 1.]
ring_locations = [0., 1.]

#=---------------------- Automated problem setup -----------------------------=#
num_rings = size(ring_particles)[1]
particles = Vorticity3DSimpleCollector()
for i = 1 : num_rings
    c = Vector3D(ring_locations[i], 0., 0.)
    n = Vector3D(1., 0., 0.)
    r = ring_radii[i]
    strength = ring_strength[i]
    n_ring_particles = ring_particles[i]
    ring = vortex_particle_ring(
        c, n, r, strength, n_ring_particles, threed_winckelmans_kernels()
        )
    particles = Vorticity3DSimpleCollector(
        get_children(particles),
        get_children(ring))
end
num_particles = size(particles)[1]

#=-------------------------- ODE time integration ----------------------------=#
for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        points = zeros(3, 0)
        point_vorticity = zeros(3, num_particles)
        cells = Array{WriteVTK.MeshCell, 1}(undef, 0)
        for j = 1 : num_particles
            points, cells = add_to_VtkMesh(points, cells, particles[j].geometry)
            point_vorticity[:, j] = [particles[j].vorticity.x,
                particles[j].vorticity.y, particles[j].vorticity.z]
        end
        vtkfile = WriteVTK.vtk_grid(string(basepath, i), points, cells)
        WriteVTK.vtk_cell_data(vtkfile, point_vorticity, "vorticity")
        outfiles = WriteVTK.vtk_save(vtkfile)
    end

    # Calculate the next iteration
    euler!(particles, particles, dt)
end
#=------------------- Now repeat until it doesn't blow up --------------------=#
end

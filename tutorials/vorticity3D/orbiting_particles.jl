#===============================================================================
orbiting_particles.jl

2 vortex particles orbiting each other - as simple as it gets.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow.
include("../../src/lowOrder3D/Vorticity3DSimpleCollector.jl")
include("../../src/lowOrder3D/VortexParticle3D.jl")
import WriteVTK  # We'll use this package to output to VTK for visualisation.

#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
const num_steps = 2000
const dt = 0.1
# Data saving parameters
basepath = "./output/orbiting_particles_"   # Where to write our output files
const save_every = 10                    # Save every 10 steps.
# Initial conditions
particle1 = VortexParticle3D(
    Vector3D(-1., 0., 0.),  # Coordinates
    Vector3D(0., 0., 1.),   # Vorticity
    0.1,                        # Radius
    threed_gaussian_kernels()
)

particle2 = deepcopy(particle1)
particle2.coord = Vector3D([1., 0., 0.])      # We can also assign to our second particle.
particle2.vorticity = convert(Vector3D, [0., 0., 0.5])
particle2.radius = 0.08

#=-------------------------- ODE time integration ----------------------------=#
particles = Vorticity3DSimpleCollector(particle1, particle2)
num_particles = length(particles)

for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        points = zeros(3, num_particles)
        point_vorticity = zeros(3, num_particles)
        cells = Array{WriteVTK.MeshCell, 1}(undef, num_particles)
        for j = 1 : num_particles
            points[:,j] = [particles[j].coord.x,
                particles[j].coord.y, particles[j].coord.z]
            point_vorticity[:, j] = [particles[j].vorticity.x,
                particles[j].vorticity.y, particles[j].vorticity.z]
            cells[j] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [j])
        end
        vtkfile = WriteVTK.vtk_grid(string(basepath, i), points, cells)
        WriteVTK.vtk_point_data(vtkfile, point_vorticity, "vorticity")
        outfiles = WriteVTK.vtk_save(vtkfile)
    end

    # Calculate the next iteration
    euler!(particles, particles, dt)
end

#=------------------- Now repeat until it doesn't blow up --------------------=#

#===============================================================================
orbiting_particles.jl

2 vortex particles orbiting each other - as simple as it gets.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow.
include("../../src/lowOrder3D/ThreeDVorticitySimpleCollector.jl")
include("../../src/lowOrder3D/ThreeDVortexParticle.jl")
using WriteVTK  # We'll use this package to output to VTK for visualisation.

#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
const num_steps = 2000
const dt = 0.1
# Data saving parameters
basepath = "./output/orbiting_particles_"   # Where to write our output files
const save_every = 10                    # Save every 10 steps.
# Initial conditions
particle1 = ThreeDVortexParticle(
    ThreeDVector(-1., 0., 0.),  # Coordinates
    ThreeDVector(0., 0., 1.),   # Vorticity
    0.1,                        # Radius
    threed_gaussian_kernels()
)

particle2 = deepcopy(particle1)
particle2.coord = ThreeDVector([1., 0., 0.])      # We can also assign to our second particle.
particle2.vorticity = convert(ThreeDVector, [0., 0., 0.5])
particle2.radius = 0.08

#=-------------------------- ODE time integration ----------------------------=#
particles = ThreeDVorticitySimpleCollector(particle1, particle2)
num_particles = length(particles)

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

#=------------------- Now repeat until it doesn't blow up --------------------=#

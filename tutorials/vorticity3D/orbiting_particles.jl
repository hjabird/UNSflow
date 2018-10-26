#===============================================================================
orbiting_particles.jl

2 vortex particles orbiting each other - as simple as it gets.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
push!(LOAD_PATH, "../../src")
import UNSflow

let
#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
num_steps = 2000
dt = 0.1
# Data saving parameters
basepath = "./output/orbiting_particles_"   # Where to write our output files
save_every = 10                    # Save every 10 steps.
# Initial conditions
particle1 = UNSflow.VortexParticle3D(
    UNSflow.Vector3D(-1., 0., 0.),  # Coordinates
    UNSflow.Vector3D(0., 0., 1.),   # Vorticity
    0.1,                        # Radius
    UNSflow.threed_gaussian_kernels()
)

particle2 = deepcopy(particle1)
# We can also assign to our second particle.
particle2.geometry.coord = UNSflow.Vector3D([1., 0., 0.])     
particle2.vorticity = convert(UNSflow.Vector3D, [0., 0., 0.5])
particle2.radius = 0.08

#=-------------------------- ODE time integration ----------------------------=#
particles = UNSflow.Vorticity3DSimpleCollector(particle1, particle2)
num_particles = length(particles)

for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        mesh = UNSflow.UnstructuredMesh()
        push!(mesh, particles)
        UNSflow.to_vtk_file(mesh, string(basepath, i))
    end

    # Calculate the next iteration
    UNSflow.euler!(particles, particles, dt)
end

#=------------------- Now repeat until it doesn't blow up --------------------=#
end #let

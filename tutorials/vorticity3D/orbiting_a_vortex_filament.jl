#===============================================================================
orbiting_a_vortex_filament.jl

A few particles orbiting round a vortex filament.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
push!(LOAD_PATH, "../../src/")
import UNSflow
include("VortexFlowFeatures.jl")

let
#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
num_steps = 2000
dt = 0.01
# Data saving parameters
basepath = "./output/orbiting_a_filament_"   # Where to write our output files
save_every = 10                       # Save every 10 steps.
# Initial conditions
filament = UNSflow.StraightVortexFilament(
    UNSflow.Vector3D(0., 0., 0.), # Start
    UNSflow.Vector3D(1., 0., 0.), # End
    0.5 # Strength
)
particle1 = UNSflow.VortexParticle3D(
    UNSflow.Vector3D(0.5, 0.5, 0.),  # Coordinates
    UNSflow.Vector3D(0., 0., 0.3),   # Vorticity
    0.1                        # Radius
)

particle2 = deepcopy(particle1)
# We can also assign to our second particle.
particle2.geometry.coord = convert(UNSflow.Vector3D, [0.5, 0., 0.5])
particle2.vorticity = convert(UNSflow.Vector3D, [0., 0.3, 0.0])
particle2.radius = 0.1

#=-------------------------- ODE time integration ----------------------------=#
particles = UNSflow.Vorticity3DSimpleCollector(particle1, particle2)
all_bodies = UNSflow.Vorticity3DSimpleCollector(filament, particles)

num_particles = length(particles)
for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        mesh = UNSflow.UnstructuredMesh()
        push!(mesh, all_bodies)
        UNSflow.to_vtk_file(mesh, string(basepath, i))
    end

    # Calculate the next iteration.
    # We are convecting particles due to the effects of all bodies.
    UNSflow.euler!(particles, all_bodies, dt)
end
end #let
#=------------------- Now repeat until it doesn't blow up --------------------=#

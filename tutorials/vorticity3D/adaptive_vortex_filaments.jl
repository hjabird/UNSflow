#===============================================================================
adaptive_vortex_filaments.jl

Showing Crow instability using  some adaptive vortex filaments.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
include("../../src/lowOrder3D/Vorticity3DSimpleCollector.jl")
include("../../src/lowOrder3D/VortexParticle3D.jl")
include("../../src/lowOrder3D/VortexParticleFilamentAdaptive.jl")
include("../../src/lowOrder3D/DiscreteGeometry3DToVTK.jl")
include("VortexFlowFeatures.jl")
import WriteVTK  # We'll use this package to output to VTK for visualisation.
import ForwardDiff
#-------------------------- User (that's you) parameters ---------------------=#
let
# ODE integration parameters
num_steps = 200
dt = 0.01
# Data saving parameters
basepath = "./output/crow_instability_"   # Where to write our output files
save_every = 1                     # Save every 10 steps
# Adaptive update parameter
update_adaptivity_every = 5        # Update our representation every 10
# Initial conditions
filament1 = VortexParticleFilamentAdaptive(
    x::Float64 -> Vector3D(0, 0, 3*x), 0.2, 0.05)
filament2 =VortexParticleFilamentAdaptive(
    x::Float64 -> Vector3D(0.5 + 0.01 * sin(pi * x), 0, 3*x), -0.2, 0.05)
collector = Vorticity3DSimpleCollector(filament1, filament2)

#=-------------------------- ODE time integration ----------------------------=#
for i = 1 : num_steps
    num_particles = length(get_children_recursive(collector))
    println("Step ", i, " with ", num_particles, " particles. ")
    if (i - 1) % update_adaptivity_every == 0
        for child in get_children(collector, Vorticity3DAdaptive)
            println("precall")
            adaptive_update!(child)
            println("postcall")
        end
    end

    # Save the current state to vtk if required
    num_particles = length(get_children_recursive(collector))
    if (i - 1) % save_every == 0
        points = zeros(3, 0)
        point_vorticity = zeros(3, num_particles)
        cells = Vector{WriteVTK.MeshCell}(undef, 0)
        particles = get_children_recursive(collector)
        for j = 1 : num_particles
            points, cells = add_to_VtkMesh(points, cells, particles[j].geometry)
            point_vorticity[:, j] =
            [   particles[j].vorticity.x,
                particles[j].vorticity.y,
                particles[j].vorticity.z    ]
        end
        vtkfile = WriteVTK.vtk_grid(string(basepath, i), points, cells)
        WriteVTK.vtk_point_data(vtkfile, point_vorticity, "vorticity")
        outfiles = WriteVTK.vtk_save(vtkfile)
    end

    # Calculate the next iteration.
    # We are convecting particles due to the effects of all bodies.
    euler!(collector, collector, dt)
end
end
#=------------------- Now repeat until it doesn't blow up --------------------=#

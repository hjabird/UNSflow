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
include("VortexFlowFeatures.jl")
import WriteVTK  # We'll use this package to output to VTK for visualisation.
import ForwardDiff
#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
const num_steps = 2000
const dt = 0.01
# Data saving parameters
basepath = "./output/crow_instability_"   # Where to write our output files
const save_every = 10                     # Save every 10 steps
# Adaptive update parameter
const update_adaptivity_every = 10        # Update our representation every 10
# Initial conditions
filament1 = VortexParticleFilamentAdaptive(
    x::Float64 -> Vector3D(0, 0, x), 1.0, 0.05)
filament2 =VortexParticleFilamentAdaptive(
    x::Float64 -> Vector3D(0.5 + 0.01 * sin(4 * pi * x), 0, x), -1.0, 0.05)
collector = Vorticity3DSimpleCollector(filament1, filament2)

#=-------------------------- ODE time integration ----------------------------=#
num_particles = length(particles)
for i = 1 : num_steps
    if (i - 1) % update_adaptivity_every == 0
        map(x->adaptive_update!(x), get_children(collector, Vorticity3DAdaptive))
    end

    # Save the current state to vtk if required
    num_particles = length(particles)
    if (i - 1) % save_every == 0
        points = zeros(3, num_particles + 2)
        point_vorticity = zeros(3, num_particles + 2)
        cells = Vector{WriteVTK.MeshCell}(undef, num_particles + 1)
        particles = get_children_recursive(collector)
        for j = 1 : num_particles
            points[:,j] =
            [   particles.children[j].coord.x,
                particles.children[j].coord.y,
                particles.children[j].coord.z    ]
            point_vorticity[:, j] =
            [   particles.children[j].vorticity.x,
                particles.children[j].vorticity.y,
                particles.children[j].vorticity.z    ]
            cells[j] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [j])
        end
        points[:, num_particles + 1] = convert(Array{Float64,1}, filament.start_coord)
        points[:, num_particles + 2] = convert(Array{Float64,1}, filament.end_coord)
        point_vorticity[:, num_particles + 1] = convert(Array{Float64,1},
            unit(filament.end_coord - filament.start_coord) * filament.vorticity)
        point_vorticity[:, num_particles + 2] =
            point_vorticity[:, num_particles + 1]
        cells[num_particles + 1] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE,
            [num_particles + 1, num_particles + 2])
        vtkfile = WriteVTK.vtk_grid(string(basepath, i), points, cells)
        WriteVTK.vtk_point_data(vtkfile, point_vorticity, "vorticity")
        outfiles = WriteVTK.vtk_save(vtkfile)
    end

    # Calculate the next iteration.
    # We are convecting particles due to the effects of all bodies.
    euler!(particles, all_bodies, dt)
end

#=------------------- Now repeat until it doesn't blow up --------------------=#

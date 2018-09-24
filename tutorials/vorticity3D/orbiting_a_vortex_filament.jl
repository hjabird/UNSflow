#===============================================================================
orbiting_a_vortex_filament.jl

A few particles orbiting round a vortex filament.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
include("../../src/lowOrder3D/Vorticity3DSimpleCollector.jl")
include("../../src/lowOrder3D/VortexParticle3D.jl")
include("../../src/lowOrder3D/StraightVortexFilament.jl")
include("ThreeDVortexParticleFlowFeatures.jl")
import WriteVTK  # We'll use this package to output to VTK for visualisation.

#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
const num_steps = 2000
const dt = 0.01
# Data saving parameters
basepath = "./output/orbiting_a_filament_"   # Where to write our output files
const save_every = 10                       # Save every 10 steps.
# Initial conditions
filament = StraightVortexFilament(
    Vector3D(0., 0., 0.), # Start
    Vector3D(1., 0., 0.), # End
    0.5 # Strength
)
particle1 = VortexParticle3D(
    Vector3D(0.5, 0.5, 0.),  # Coordinates
    Vector3D(0., 0., 0.3),   # Vorticity
    0.1                        # Radius
)

particle2 = deepcopy(particle1)
# We can also assign to our second particle.
particle2.coord = convert(Vector3D, [0.5, 0., 0.5])
particle2.vorticity = convert(Vector3D, [0., 0.3, 0.0])
particle2.radius = 0.1

#=-------------------------- ODE time integration ----------------------------=#
particles = Vorticity3DSimpleCollector(particle1, particle2)
all_bodies = Vorticity3DSimpleCollector(filament, particles)

num_particles = length(particles)
for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        points = zeros(3, num_particles + 2)
        point_vorticity = zeros(3, num_particles + 2)
        cells = Vector{WriteVTK.MeshCell}(undef, num_particles + 1)
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

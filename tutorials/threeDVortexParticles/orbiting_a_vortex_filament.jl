#===============================================================================
orbiting_a_vortex_filament.jl

A few particles orbiting round a vortex filament.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow.
using Dierckx
using SpecialFunctions
import Base.size
include("../../src/kinem.jl")
include("../../src/utils.jl")
include("../../src/delVort.jl")
include("../../src/lowOrder2D/typedefs.jl")
include("../../src/lowOrder2D/calcs.jl")
include("../../src/lowOrder3D/typedefs.jl")
include("../../src/lowOrder3D/calcs.jl")
using WriteVTK  # We'll use this package to output to VTK for visualisation.

#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
const num_steps = 2000
const dt = 0.01
# Options: euler_forward_step!, explicit_midpoint_step!
ode_method = explicit_midpoint_step!
# Data saving parameters
basepath = "./output/orbiting_a_filament_"   # Where to write our output files
const save_every = 10                       # Save every 10 steps.
# Initial conditions
filament = ThreeDStraightVortexFilament(
    ThreeDVector(0., 0., 0.), # Start
    ThreeDVector(1., 0., 0.), # End
    0.5 # Strength
)
particle1 = ThreeDVortexParticle(
    ThreeDVector(0.5, 0.5, 0.),  # Coordinates
    ThreeDVector(0., 0., 0.3),   # Vorticity
    0.1,                        # Radius
    # These don't matter, but have to be filled in anyway.
    ThreeDVector(-1., 0., 0.),  # Velocity
    ThreeDVector(0., 0., 1.)    # Rate of change of vorticity
)

particle2 = deepcopy(particle1)
# We can also assign to our second particle.
particle2.coord = convert(ThreeDVector, [0.5, 0., 0.5])
particle2.vorticity = convert(ThreeDVector, [0., 0.3, 0.0])
particle2.size = 0.1

# options: singular, planetary, exponential, winckelmans, tanh, gaussian,
# and super_gaussian
reduction_factor_fn, vorticity_fraction_fn = threed_super_gaussian_kernels()

#=-------------------------- ODE time integration ----------------------------=#
particle_set = ThreeDVortexParticleSet(
    [particle1, particle2],
    reduction_factor_fn,
    vorticity_fraction_fn
)
function fil_ind_vel(particle :: ThreeDVortexParticle)
    return ind_vel(filament, particle.coord)
end
function fil_ind_dvortdt(particle :: ThreeDVortexParticle)
    return ind_dvortdt(particle, filament)
end

num_particles = size(particle_set.particles)[1]
for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        points = zeros(3, num_particles + 2)
        point_vorticity = zeros(3, num_particles + 2)
        cells = Array{MeshCell, 1}(num_particles + 1)
        for j = 1 : num_particles
            points[:,j] =
            [   particle_set.particles[j].coord.x,
                particle_set.particles[j].coord.y,
                particle_set.particles[j].coord.z    ]
            point_vorticity[:, j] =
            [   particle_set.particles[j].vorticity.x,
                particle_set.particles[j].vorticity.y,
                particle_set.particles[j].vorticity.z    ]
            cells[j] = MeshCell(VTKCellTypes.VTK_VERTEX, [j])
        end
        points[:, num_particles + 1] = convert(Array{Float64,1}, filament.start_coord)
        points[:, num_particles + 2] = convert(Array{Float64,1}, filament.end_coord)
        point_vorticity[:, num_particles + 1] = convert(Array{Float64,1},
            unit(filament.end_coord - filament.start_coord) * filament.vorticity)
        point_vorticity[:, num_particles + 2] =
            point_vorticity[:, num_particles + 1]
        cells[num_particles + 1] = MeshCell(VTKCellTypes.VTK_LINE,
            [num_particles + 1, num_particles + 2])
        vtkfile = vtk_grid(string(basepath, i), points, cells)
        vtk_point_data(vtkfile, point_vorticity, "vorticity")
        outfiles = vtk_save(vtkfile)
    end

    # Calculate the next iteration
    particles = ode_method(particle_set, fil_ind_vel, fil_ind_dvortdt, dt)
end

#=------------------- Now repeat until it doesn't blow up --------------------=#

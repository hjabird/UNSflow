#===============================================================================
state_space_sheet_rollup.jl

A demo of a sheet of vortex particles rollling up, but this time using
the DifferenetialEquations package and a state space representation to get
the system solved more accurately and, hopefully, more quickly.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow instead.
let
include("../../src/lowOrder3D/VortexParticle3D.jl")
include("../../src/lowOrder3D/Vorticity3D.jl")
include("../../src/lowOrder3D/Vorticity3DSimpleCollector.jl")
include("../../src/lowOrder3D/Vortex3DRegularisationFunctions.jl")
include("../../src/lowOrder3D/DiscreteGeometry3DToVTK.jl")
include("VortexFlowFeatures.jl")
import WriteVTK  # We'll use this package to output to VTK for visualisation.
import DifferentialEquations

#---------------------------- User parameters --------------------------------=#
# ODE integration parameters
max_time = 15                         # seconds.
# Data saving parameters
basepath = "./output/state_space_sheet_rollup_"   # Where to write our output files
save_every = 0.5                      # Save every 0.5 seconds

# We can make a sheet of vortex particles defined by functions. Lets make
# a flat sheet for now, because the function that does the work here
# currently only works for nearly flat stuff.
# Define a function that turns x,y to TheeDVector
function f1(x, y)
    return Vector3D(x, y, 0)
end
# We also want it to be bounded
bounds = [-1, 2, -0.5, 0.5] # minx, maxx, miny, maxy
# And it needs to have a known number of particles
np_x = 40
np_y = 15
# and finally we define a continious vorticy density function function.
# In theory, out of plane vorticity is not ok unless we have thickness, but
# we'll conveniently gloss over that. We define vorticity here in global coord.
function f2(x, y)
    return Vector3D(0.05 * y, 0.1 * (x^2)^0.25, 0.0)
end

#=---------------------- Automated problem setup -----------------------------=#
particles = vortex_particle_sheet(
    f1,
    f2,
    bounds,
    np_x,
    np_y,
    threed_winckelmans_kernels()
)
num_particles = length(particles)

#=-------------------------- ODE time integration ----------------------------=#
time = 0.0
state = state_vector(particles)
state_deriv_calls = 0
function dstate(state_vect :: Vector{Float64})
    update_using_state_vector!(particles, state_vect)
    state_deriv_calls += 1
    return state_time_derivative(particles, particles)
end
f(u, p, t) = dstate(u)
while time < max_time
    update_using_state_vector!(particles, state)
    points = zeros(3, 0)
    point_vorticity = zeros(3, num_particles)
    cells = Array{WriteVTK.MeshCell, 1}(undef, 0)
    for j = 1 : num_particles
        points, cells = add_to_VtkMesh(points, cells, particles[j].geometry)
        point_vorticity[:, j] = [particles[j].vorticity.x,
            particles[j].vorticity.y, particles[j].vorticity.z]
    end
    vtkfile = WriteVTK.vtk_grid(string(basepath, Int32(time/save_every)),
        points, cells)
    WriteVTK.vtk_point_data(vtkfile, point_vorticity, "vorticity")
    outfiles = WriteVTK.vtk_save(vtkfile)


    # Calculate the next iteration
    prob = DifferentialEquations.ODEProblem(f, state, (time, time+save_every))
    solution = DifferentialEquations.solve(prob)
    state = solution.u[end]
    time += save_every
end
println("Called ind_vel ", state_deriv_calls, " times.")
#=------------------- Now repeat until it doesn't blow up --------------------=#
end

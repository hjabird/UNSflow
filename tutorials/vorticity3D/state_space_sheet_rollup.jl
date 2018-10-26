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

import DifferentialEquations
push!(LOAD_PATH,"../../src/")
import UNSflow
include("VortexFlowFeatures.jl")
let
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
    return UNSflow.Vector3D(x, y, 0)
end
# We also want it to be bounded
bounds = [-1, 2, -0.5, 0.5] # minx, maxx, miny, maxy
# And it needs to have a known number of particles
np_x = 30
np_y = 10
# and finally we define a continious vorticy density function function.
# In theory, out of plane vorticity is not ok unless we have thickness, but
# we'll conveniently gloss over that. We define vorticity here in global coord.
function f2(x, y)
    return UNSflow.Vector3D(0.05 * y, 0.1 * (x^2)^0.25, 0.0)
end

#=---------------------- Automated problem setup -----------------------------=#
particles = vortex_particle_sheet(
    f1,
    f2,
    bounds,
    np_x,
    np_y,
    UNSflow.threed_winckelmans_kernels()
)
num_particles = length(particles)

#=-------------------------- ODE time integration ----------------------------=#
time = 0.0
state = UNSflow.state_vector(particles)
state_deriv_calls = 0
function dstate(state_vect :: Vector{Float64})
    UNSflow.update_using_state_vector!(particles, state_vect)
    state_deriv_calls += 1
    return UNSflow.state_time_derivative(particles, particles)
end
f(u, p, t) = dstate(u)
while time < max_time
    UNSflow.update_using_state_vector!(particles, state)
    mesh = UNSflow.UnstructuredMesh()
    push!(mesh, particles)
    UNSflow.to_vtk_file(mesh, string("output/state_space_sheet_rollup_", 
        Int64(round(time / save_every))))

    # Calculate the next iteration
    prob = DifferentialEquations.ODEProblem(f, state, (time, time+save_every))
    @time solution = DifferentialEquations.solve(prob)
    state = solution.u[end]
    time += save_every
end
println("Called ind_vel ", state_deriv_calls, " times.")
#=------------------- Now repeat until it doesn't blow up --------------------=#
end

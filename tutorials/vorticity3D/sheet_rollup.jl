#===============================================================================
sheet_rollup.jl

A demo of a sheet of vortex particles rollling up.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
push!(LOAD_PATH,"../../src/")
import UNSflow
include("VortexFlowFeatures.jl")

let
#---------------------------- User parameters --------------------------------=#
# ODE integration parameters
num_steps = 300
dt = 0.05
# Data saving parameters
basepath = "./output/sheet_rollup_"   # Where to write our output files
save_every = 10                       # Save every 10 steps.

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
np_x = 40
np_y = 15
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
for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        mesh = UNSflow.UnstructuredMesh()
        push!(mesh, particles)
        UNSflow.to_vtk_file(mesh, string("output/sheet_rollup_", i))
    end

    # Calculate the next iteration
    UNSflow.euler!(particles, particles, dt)
end

#=------------------- Now repeat until it doesn't blow up --------------------=#
end #let

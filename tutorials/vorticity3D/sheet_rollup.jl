#===============================================================================
sheet_rollup.jl

A demo of a sheet of vortex particles rollling up.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow instead.
let
import WriteVTK  # We'll use this package to output to VTK for visualisation.
push!(LOAD_PATH,"../../src/")
import UNSflow
include("VortexFlowFeatures.jl")


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
        points = zeros(3, 0)
        point_vorticity = zeros(3, num_particles)
        cells = Array{WriteVTK.MeshCell, 1}(undef, 0)
        for j = 1 : num_particles
            points, cells = UNSflow.add_to_VtkMesh(points, cells, particles[j].geometry)
            point_vorticity[:, j] = [particles[j].vorticity.x,
                particles[j].vorticity.y, particles[j].vorticity.z]
        end
        vtkfile = WriteVTK.vtk_grid(string(basepath, i), points, cells)
        WriteVTK.vtk_point_data(vtkfile, point_vorticity, "vorticity")
        outfiles = WriteVTK.vtk_save(vtkfile)
    end

    # Calculate the next iteration
    UNSflow.euler!(particles, particles, dt)
end

#=------------------- Now repeat until it doesn't blow up --------------------=#
end

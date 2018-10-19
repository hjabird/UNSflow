#===============================================================================
    steady_vortex_lattice_method.jl

    An example of using the vortex lattice method.

    Initial code: HJAB 2018

    Copyright (c) 2018 HJA Bird

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to
    deal in the Software without restriction, including without limitation the
    rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
    sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
    IN THE SOFTWARE.
------------------------------------------------------------------------------=#
push!(LOAD_PATH, "../../src/")
import UNSflow
import WriteVTK

let
# STEP 1 : Define a wing. ------------------------------------------------------
# This is done by using equation surf. We take x=-1 as the leading edge, x=1
# as the trailing edge. 
aoa = 0.1
aspect_ratio = 4
# Rectangular:
#surf_fn = x->UNSflow.Vector3D(x[1]/2, x[2]*aspect_ratio/2, 0)
# Or perhaps elliptic:
surf_fn = x->UNSflow.Vector3D(sqrt(1- x[2]^2) * x[1]/2, x[2]*aspect_ratio/2, 0)
# Apply the angle of attack:
surf_fn_aoa = x->UNSflow.rotate_about_y(surf_fn(x), aoa)
# The generate the surface and discretise the surface:
wing_surf = UNSflow.EquationSurf(surf_fn_aoa)
spanwise_disc = map(x->cbrt(x), -1:0.1:1)
chordwise_disc = map(x->x, -1 : 0.125 : 1)
wing_geom = UNSflow.discretise(wing_surf, UNSflow.BilinearQuadSurf, 
                                                chordwise_disc, spanwise_disc)

# STEP 2: Define the free stream and wake parameters to setup the problem ------
free_stream = UNSflow.Vector3D(1., 0, 0)
# We have to decide on how to discretise the wake as it travels downstream.
downstream = map(x->x^1.4, collect(0:0.1:3)[2:end] .* 5)
problem = UNSflow.VortexLatticeMethod(wing_geom, free_stream, downstream)

# STEP 3: Solve the system -----------------------------------------------------
for i = 1 : 1
    # We call solve to calculate vorticities
    UNSflow.solve!(problem)
    # And then relax_wake! to let the wake curl up as expected.
    UNSflow.relax_wake!(problem; relaxation_factor=0.5)
    # We might want to solve again now the problem geometry has changed...
end

# STEP 4: Post-process ---------------------------------------------------------
# We need the induced velocity to caluclate the forces on the lattice:
ind_vel_fn = x -> free_stream + 
    UNSflow.induced_velocity(problem.wake_aerodynamic, x) + 
    UNSflow.induced_velocity(problem.wing_aerodynamic, x)
# We can calculate some useful stuff using steady_loads:
force, moment, pressure = UNSflow.steady_loads(
    problem.wing_aerodynamic, 
    ind_vel_fn; 
    # The wing vortex lattice is overlapped by that of the wake, so we need to 
    # correct for that - imax indicates the maximum i index of the lattice.
    imax_filament_strs=-problem.wake_aerodynamic.strengths[1,:])

# Print out some useful things:
println("Wing area is ", UNSflow.area(wing_geom))
println("Forces are ", force)
println("Force coeffs ", 2 * force / 
                        (abs(free_stream)^2 * UNSflow.area(wing_geom)))
println("Moments are ", moment)

# STEP 5: Output to VTK file ---------------------------------------------------
pressure = vcat(vec(pressure), 
                zeros(length(problem.wake_aerodynamic.strengths)))

# And plot the shape of stuff:
points, cells = UNSflow.to_VtkMesh(
    convert(Vector{UNSflow.BilinearQuad}, problem.wing_geometry))
points, cells = UNSflow.add_to_VtkMesh(
    points, cells, 
    convert(Vector{UNSflow.BilinearQuad}, problem.wake_aerodynamic.geometry))
vtkfile = WriteVTK.vtk_grid("output/steady_vortex_lattice", points, cells)
vorticity = vcat(
    vec(problem.wing_aerodynamic.strengths), 
    vec(problem.wake_aerodynamic.strengths))
WriteVTK.vtk_cell_data(vtkfile, vorticity, "vorticity")
WriteVTK.vtk_cell_data(vtkfile, pressure, "pressure")
outfiles = WriteVTK.vtk_save(vtkfile)

end #let

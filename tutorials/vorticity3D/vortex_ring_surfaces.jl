#===============================================================================
    vortex_ring_surfaces.jl

    An example of using a vortex lattice.

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

import WriteVTK
push!(LOAD_PATH,"../../src/")
import UNSflow

let 
# This is an excuse to use a vortex lattice surface in as simple a context
# as possible.

# Setting up a problem ---------------------------------------------------------
# First, we'll define a square vortex lattice...
vort_surf_eqn = UNSflow.EquationSurf(x->UNSflow.Vector3D(x[1], x[2], 0))
xdiv = collect(-1:2/10:1)
ydiv = collect(-1:2/10:1)
disc_surf = UNSflow.discretise(
    vort_surf_eqn, 
    UNSflow.BilinearQuadSurf, xdiv, ydiv)
# And we can make a surface (with zero vorticity) using that geometry:
ring_lattice = UNSflow.VortexRingLattice(
    zeros(length(xdiv)-1, length(ydiv)-1), # Preset vortex ring strengths
    disc_surf)  # The surface we want to base the geometry off.

# We want it to do something... So lets put it inside a vortex ring:
ext_ring = UNSflow.VortexRing(
    UNSflow.Vector3D(-2,-2,0), UNSflow.Vector3D(2,-2,0), 
    UNSflow.Vector3D(2,2,0), UNSflow.Vector3D(-2,2,0), 
    1.0 # The strength
)
# Perhaps also a good old free stream too - Enjoy D'Alembert's paradox:
free_stream = UNSflow.Vector3D(0, 0, 1)

# Solving it -------------------------------------------------------------------
# And apply some Neumann BCs:
surf_coords = vec(UNSflow.centres(disc_surf))
surf_normals = vec(UNSflow.normals(disc_surf))
mat = UNSflow.normal_velocity_influence_matrix(
    ring_lattice, vec(surf_coords), vec(surf_normals))
known = map(
    x->UNSflow.dot(UNSflow.induced_velocity(ext_ring, x[1]) +free_stream, x[2]), 
    zip(surf_coords, surf_normals))
# And solve for them:
solution = mat \ -known
UNSflow.update_using_vorticity_vector!(ring_lattice, solution)

# Doing a little post-processing -----------------------------------------------
ind_vel_fn = x -> UNSflow.induced_velocity(ring_lattice, x) +
                  UNSflow.induced_velocity(ext_ring, x) + free_stream
# We'll use density == 1
# Keyword arg measurement_centre defines where to measure moments about.
force, moment, pressure = UNSflow.steady_loads(ring_lattice, ind_vel_fn, 1.;
                                measurement_centre=UNSflow.Vector3D(0, 0, 0))
println("Plate area was ", UNSflow.area(disc_surf))
println("Moments were ", moment)

# And outputting the data to VTK -----------------------------------------------
points, cells = UNSflow.to_VtkMesh(
    convert(Vector{UNSflow.BilinearQuad}, disc_surf))
points, cells = UNSflow.add_to_VtkMesh(points, cells, ext_ring.geometry)
# Get the celldata into the required format.
vorticity = vec(ring_lattice.strengths)
push!(vorticity, 1.0) #For our vortex ring.
pressure = vec(pressure)
push!(pressure, 0.0) #For our vortex ring

vtkfile = WriteVTK.vtk_grid("output/lattice_in_ring", points, cells)
WriteVTK.vtk_cell_data(vtkfile, vorticity, "vorticity")
WriteVTK.vtk_cell_data(vtkfile, pressure, "pressure")
outfiles = WriteVTK.vtk_save(vtkfile)

end # Let
# That's all folks -------------------------------------------------------------

#===============================================================================
    G4_using_surface_properties.jl

    We'll use EquationSurf.derivative and EquationSurf.normal to
    offset our surface.

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

push!(LOAD_PATH,"../../src/")
import UNSflow

let
#= Code to make a tube stolen from G1_making_a_tube.jl ========================#
z_def = x->x[2]
x_def = x->sin(pi * x[1])
y_def = x->cos(pi * x[1])
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
#= Finish making a tube =======================================================#

# MODIFYING A SURFACE USING ITS NORMAL VECTOR
# We want a surface offset from our original one by a variable amount.
# Lets define this offset as an equation:
offset = x -> abs(0.75 * sin(pi * 2 * x[1]) * cos(pi * x[2] / 2))
# Now we generate an equation for the surface position...
xyz_def = x -> UNSflow.evaluate(surf, x) - UNSflow.normal(surf, x) * offset(x)
# And use that to generate our new surface:
offset_surf = UNSflow.EquationSurf(xyz_def)
disc_range = collect(-1:0.02:1)
discrete_surf = UNSflow.discretise(offset_surf, UNSflow.BilinearQuad,
    disc_range, disc_range)

# EXTRACTING THE TANGENTS FROM A SURFACE

# We might also want to know, for example, the tangents of our surface.
# First, we want the centre of each of the quads we discretised into:
tang_range = (disc_range[1:end-1] .+ disc_range[2:end]) / 2
# Its worth knowing here that Eqsurf discretises in the same way that a Matrix
# is converted to linear indexing.
tangent_coords = vec(reshape([[x, y] for x=tang_range, y=tang_range], :, 1))
# And generate our tangents:
td1 = map(x->UNSflow.derivative(offset_surf, 1,x), tangent_coords)
td2 = map(x->UNSflow.derivative(offset_surf, 2,x), tangent_coords)
# We can also extract element areas:
areas = vec(UNSflow.area.(discrete_surf))

# And save to VTK
mesh = UNSflow.UnstructuredMesh()
# Add cells returns a vector of [(cell_idx::Int64, cell_vertex_idxs[])]
percell_cidx_pidx = UNSflow.add_cells!(mesh, discrete_surf)
map(x->UNSflow.add_celldata!(mesh, x[1][1], "Tangentd1", x[2]),
    zip(percell_cidx_pidx, td1))
map(x->UNSflow.add_celldata!(mesh, x[1][1], "Tangentd2", x[2]),
    zip(percell_cidx_pidx, td2))
map(x->UNSflow.add_celldata!(mesh, x[1][1], "Area", x[2]),
    zip(percell_cidx_pidx, areas))
UNSflow.to_vtk_file(mesh, "output/G4_using_surface_props")

end #let

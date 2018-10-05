#===============================================================================
    equation_surf_surf_properties.jl

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
import WriteVTK

let
#= Code to make a tube stolen from equation_surf_tube.jl ======================#
z_def = x->x[2]
x_def = x->sin(pi * x[1])
y_def = x->cos(pi * x[1])
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
#= Finish making a tube =======================================================#

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
# We might also want to know, for example, the tangents of our surface.
# First, we want the centre of each of the quads we discretised into:
tang_range = (disc_range[1:end-1] .+ disc_range[2:end]) / 2
# Its worth knowing here that Eqsurf discretises such that the y index is
# changing the most ie: [1,1],[1,2],[1,3],[2,1],[2,2],[2,3]...
tangent_coords = vec(reshape([[x, y] for y=tang_range,x=tang_range], :, 1))
# And generate our tangents:
td1 = map(x->UNSflow.derivative(offset_surf, 1,x), tangent_coords)
td2 = map(x->UNSflow.derivative(offset_surf, 2,x), tangent_coords)
# And save to VTK
# We need to convert our vector of tangents to matrices...
td1_mat = convert(Matrix{Float64}, td1)
td2_mat = convert(Matrix{Float64}, td2)

points, cells = UNSflow.to_VtkMesh(discrete_surf)
vtkfile = WriteVTK.vtk_grid("output/equation_surf_surf_properties", points, cells)
vtkfile = WriteVTK.vtk_cell_data(vtkfile, td1_mat, "Tangentd1")
vtkfile = WriteVTK.vtk_cell_data(vtkfile, td2_mat, "Tangentd2")
outfiles = WriteVTK.vtk_save(vtkfile)

end #let

#===============================================================================
    equation_surf_tube.jl

    An example of how one might use EquationSurf to generate a surface.
    The surface is then discretised and dumped into a VTK file for veiwing.

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

# First, we're going to make a tube.
# We'll put the linear direction in z
z_def = x->x[2]
# And the cross section in x,y, remebering we have [-1,1] to work with.
x_def = x->sin(pi * x[1])
y_def = x->cos(pi * x[1])
# And now we can make our surface.
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
# To export it we need to discretise it. We'll turn it into some
# BilinearQuad elements. We'll use a coarser discretisation in the z direction.
discrete_surf = UNSflow.discretise(surf, UNSflow.BilinearQuad,
    collect(-1:0.1:1), collect(-1:0.2:1))

# ... And now we can save it to a file:
points, cells = UNSflow.to_VtkMesh(discrete_surf)
vtkfile = WriteVTK.vtk_grid("output/equation_surf_tube", points, cells)
outfiles = WriteVTK.vtk_save(vtkfile)

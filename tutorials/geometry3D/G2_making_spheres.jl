#===============================================================================
    G2_making_spheres.jl

    We'll generate a sphere and output the discretisation to a VTK file.
    Because we are going to generate two hemisphere, this is ever so slightly
    more complicated than equation_surf_tube.jl

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
# Lets define the surface of a sphere such that our x[2] is the longitude
# and x[1] is latitude:
z_def = x->sin(x[2] * pi / 2)
x_def = x->cos(x[1] * pi) * cos(x[2] * pi/2)
y_def = x->sin(x[1] * pi) * cos(x[2] * pi/2)
# And now we can make our surface.
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
# To export it we need to discretise it. We'll turn it into some
# BilinearQuad elements. We'll use a coarser discretisation in the z direction.
discrete_surf = UNSflow.discretise(surf, UNSflow.BilinearQuad,
    collect(-1:0.1:1), collect(-1:0.05:1))

# ... And now we can save it to a file:
points, cells = UNSflow.to_VtkMesh(discrete_surf)
vtkfile = WriteVTK.vtk_grid("output/G2_bad_sphere", points, cells)
outfiles = WriteVTK.vtk_save(vtkfile)

# Take a look at the file - we have a huge number of cells at the poles, and
# fewer at the equator. A classical problem. We can do better by splitting our
# sphere into six in pi / 2 by pi / 2 arcs things.

x_def = x->sqrt(2)/2 * x[1] * sqrt(1 - x[2]^2 / 2)
y_def = x->sqrt(2)/2 * x[2] * sqrt(1 - x[1]^2 / 2)
z_def = x->sqrt(abs(1 - x_def(x)^2 - y_def(x)^2))
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
discrete_surf = UNSflow.discretise(surf, UNSflow.BilinearQuad,
    collect(-1:1/6:1), collect(-1:1/6:1))
points, cells = UNSflow.to_VtkMesh(discrete_surf)
for i = 1 : 4
    z_def = x ->sin(x[2] * pi / 4)
    x_def = x->cos((x[1]+2*i) * pi/4) * cos(x[2] * pi/4)
    y_def = x->sin((x[1]+2*i) * pi/4) * cos(x[2] * pi/4)
    surf = UNSflow.EquationSurf(x_def, y_def, z_def)
    discrete_surf = UNSflow.discretise(surf, UNSflow.BilinearQuad,
        collect(-1:1/6:1), collect(-1:1/6:1))
    points, cells = UNSflow.add_to_VtkMesh(points, cells, discrete_surf)
end
x_def = x->sqrt(2)/2 * x[1] * sqrt(1 - x[2]^2 / 2)
y_def = x->sqrt(2)/2 * x[2] * sqrt(1 - x[1]^2 / 2)
z_def = x->-sqrt(abs(1 - x_def(x)^2 - y_def(x)^2))
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
discrete_surf = UNSflow.discretise(surf, UNSflow.BilinearQuad,
    collect(-1:1/6:1), collect(-1:1/6:1))
points, cells = UNSflow.add_to_VtkMesh(points, cells, discrete_surf)

vtkfile = WriteVTK.vtk_grid("output/G2_good_sphere", points, cells)
outfiles = WriteVTK.vtk_save(vtkfile)

end #let

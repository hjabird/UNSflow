#===============================================================================
    G3_equationsurf_paths.jl

    An example of how one might pick out a path on an equation surf.

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
#= Code to make a tube stolen from G1_making_a_tube.jl ========================#
z_def = x->x[2]
x_def = x->sin(pi * x[1])
y_def = x->cos(pi * x[1])
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
discrete_surf = UNSflow.discretise(surf, UNSflow.BilinearQuad,
    collect(-1:0.1:1), collect(-1:0.2:1))
points, cells = UNSflow.to_VtkMesh(discrete_surf)
#= Finish making a tube =======================================================#

# Lets make a path on the surface. We can define some local x coordinates
# and y coordinates:
xs = collect(-0.5: 0.02 : 0.7)
ys = 0.8 * sin.(xs * pi * 4) .+ 0.1
path = UNSflow.discretise(surf, UNSflow.PolyLine2, [xs ys])
points, cells = UNSflow.add_to_VtkMesh(points, cells, path)

# We can also distribute points evenly using the Spline3D class.
xs = collect(-0.5: 0.02 : 0.7)
ys = -0.7 * sin.(xs * pi * 3 .+ 0.5)
spline = UNSflow.discretise(surf, UNSflow.Spline3D, [xs ys])
vects = UNSflow.distribute(spline, 0, spline.limits[2], 300)
points, cells = UNSflow.add_to_VtkMesh(points, cells, map(UNSflow.Point3D, vects))

# Now turn it into a vtk file...
vtkfile = WriteVTK.vtk_grid("output/G3_surface_path", points, cells)
outfiles = WriteVTK.vtk_save(vtkfile)

end #let

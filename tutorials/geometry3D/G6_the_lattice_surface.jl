#===============================================================================
    G6_the_lattice_surface.jl

    So far we've used solid surfaces. Sometimes its inappropriate to describe
    goemetry with a solid surface however. We might do better to use a lattice.
    We'll do that here.

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
# We'll steal the code to make a sphere from G2 for this.
z_def = x->sin(x[2] * pi / 2)
x_def = x->cos(x[1] * pi) * cos(x[2] * pi/2)
y_def = x->sin(x[1] * pi) * cos(x[2] * pi/2)
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
# Its very easy to define a lattice surface...
discrete_surf = UNSflow.discretise(surf, UNSflow.LatticeQuadSurf,
    collect(-1:0.1:1), collect(-1:0.05:1))

# At this point, we'd normally just export it as a VTK file, but we have a 
# problem - VTK doesn't naturally do lattices. We have some options:
# Convert to BilinearQuadSurf:
bqs = convert(UNSflow.BilinearQuadSurf, discrete_surf)
# Or a load or PolyLine2 rings:
pl2rings = convert(Vector{UNSflow.PolyLine2}, discrete_surf)
# Or a load of Line2:
l2vect = convert(Vector{UNSflow.Line2}, discrete_surf)

# This can create an annoyance whilst saving data. Not a lot can be done to fix
# this beyond diversifying file formats.

# We'll export the Line2 version because it best represents the lattice:
points, cells = UNSflow.to_VtkMesh(l2vect)
vtkfile = WriteVTK.vtk_grid("output/G6_line2_sphere", points, cells)
outfiles = WriteVTK.vtk_save(vtkfile)

# We can also do the same to make a tube. Why? Because I need to test that
# the ends appear in the conversion to line2s.
z_def = x->x[2]
x_def = x->sin(pi * x[1])
y_def = x->cos(pi * x[1])
surf = UNSflow.EquationSurf(x_def, y_def, z_def)
discrete_surf = UNSflow.discretise(surf, UNSflow.LatticeQuadSurf,
    collect(-1:0.1:1), collect(-1:0.2:1))
l2vect = convert(Vector{UNSflow.Line2}, discrete_surf)
points, cells = UNSflow.to_VtkMesh(l2vect)
vtkfile = WriteVTK.vtk_grid("output/G6_line2_tube", points, cells)
outfiles = WriteVTK.vtk_save(vtkfile)

end #let

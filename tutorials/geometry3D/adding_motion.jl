#===============================================================================
    adding_motion.jl

    We'll use UNSflow.CoordinateTransform3D to add movement to our surfaces.

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
discrete_surf = UNSflow.discretise(surf, UNSflow.BilinearQuad,
    collect(-1:0.1:1), collect(-1:0.2:1))
#= Finish making a tube =======================================================#

# We have a geometry so we can add motion. For now we'll define a twisting
# motion. Our first argument, x, is a input Vector3D, and t is the current
# time.
trans = UNSflow.CoordinateTransform3D((x,t)->
                                    UNSflow.rotate_about_z(x, sin(t) * pi/2))
# Define a dt and create a complete cycle of oscillation, outputting the results
# to a VTK file.
nsteps = 200
dt = 4*pi / nsteps
UNSflow.increment!(trans, dt)
for i = 1 : nsteps
    # We need to make a copy of the original.
    newsurf = Vector{UNSflow.BilinearQuad}()
    for subsurf in discrete_surf
        # Transform the coordinates
        newcoords = trans.(UNSflow.coords(subsurf))
        push!(newsurf, UNSflow.BilinearQuad(newcoords))
    end
    # And write to a VTK file.
    points, cells = UNSflow.to_VtkMesh(newsurf)
    str = string("output/adding_motion_", i)
    vtkfile = WriteVTK.vtk_grid(str, points, cells)
    outfiles = WriteVTK.vtk_save(vtkfile)
    # And its easy to forget to
    UNSflow.increment!(trans, dt)
end
end # let

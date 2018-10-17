#===============================================================================
    Vorticity3DToVTK.jl

    Save a Vorticity3D object to a VTK file.

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

"""
Save a Vorticity3D object to a VTK file.
"""
save_to_VTK(a::Vorticity3D)
    points, cells = UNSflow.to_VtkMesh(convert(Vector{UNSflow.BilinearQuad}, disc_surf))
    points, cells = UNSflow.add_to_VtkMesh(points, cells, ext_ring.geometry)
    vorticity = vec(ring_lattice.strengths)
    vtkfile = WriteVTK.vtk_grid("output/lattice_in_ring", points, cells)
    WriteVTK.vtk_cell_data(vtkfile, vorticity, "vorticity")
    outfiles = WriteVTK.vtk_save(vtkfile)
end

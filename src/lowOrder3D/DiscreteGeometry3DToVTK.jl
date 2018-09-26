#===============================================================================
    DiscreteGeometry3DToVTK.jl

    Add geometry to a VTK file using the WriteVTK module.

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
include("DiscreteGeometry3D.jl")
include("Line2.jl")
include("Point3D.jl")
include("PolyLine2.jl")
#try
    import WriteVTK

    function add_to_VtkMesh(
        filepoints::Matrix{Float64},
        filecells::Vector{WriteVTK.MeshCell},
        geometry::DiscreteGeometry3D
        )

        new_points = zeros(3, number_of_control_points(geometry))
        new_cell_type = vtk_cell_type(geometry)
        vec3_points = coords(geometry)
        for i = 1 : length(vec3_points)
            coord = vec3_points[i]
            new_points[:, i] = [coord.x, coord.y, coord.z]
        end
        idx_offset = size(filepoints)[2]
        filepoints = hcat(filepoints, new_points)
        newcells = deepcopy(filecells)
        push!(newcells, WriteVTK.MeshCell(new_cell_type,
            Vector{Int64}(idx_offset + 1 : idx_offset + length(vec3_points))))
        filecells = newcells
        return filepoints, filecells
    end

    function vtk_cell_type(a::DiscreteGeometry3D)
        error("VTK cell type of ", typeof(a), " has not been defined.")
    end
    function vtk_cell_type(a::Point3D)
        return WriteVTK.VTKCellTypes.VTK_VERTEX
    end
    function vtk_cell_type(a::Line2)
        return WriteVTK.VTKCellTypes.VTK_LINE
    end
    function vtk_cell_type(a::PolyLine2)
        return WriteVTK.VTKCellTypes.VTK_POLY_LINE
    end
#=catch
    function add_to_VtkMesh()
        error("Could not import WriteVTK.")
    end
end=#

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

function add_to_VtkMesh(
    filepoints :: Matrix{Float64},
    filecells :: Vector{WriteVTK.MeshCell},
    geometry :: Vector{T}
    ) where T <: DiscreteGeometry3D
    @assert(size(filepoints)[1] == 3, "expected 3 by n-points array!")
    point_count = map(number_of_control_points, geometry)
    new_points = zeros(3, sum(point_count))
    point_offset = size(filepoints)[2] + 1
    npoint_offset = 1
    new_cells = Vector{WriteVTK.MeshCell}(undef, length(point_count))
    for i = 1:length(point_count)
        geom = geometry[i]
        new_cell_type = vtk_cell_type(geom)
        coordinates = coords(geom)
        new_cells[i] = WriteVTK.MeshCell(new_cell_type,
            Vector{Int64}(point_offset : point_offset + point_count[i] - 1))
        for coordinate in coordinates
            new_points[:, npoint_offset] = Base.convert(Vector{Float64},
                    coordinate)
            npoint_offset += 1
            point_offset += 1
        end
    end
    filepoints = hcat(filepoints, new_points)
    filecells = deepcopy(filecells)
    filecells = append!(filecells, new_cells)
    return filepoints, filecells
end

function add_to_VtkMesh(
    filepoints :: Matrix{Float64},
    filecells :: Vector{WriteVTK.MeshCell},
    geometry :: Matrix{T}
    ) where T <: DiscreteGeometry3D
    return add_to_VtkMesh(filepoints, filecells, vec(geometry))
end

function to_VtkMesh(
    geometry :: Vector{T}
    ) where T <: DiscreteGeometry3D

    filepoints = zeros(3, 0)
    filecells = Vector{WriteVTK.MeshCell}(undef, 0)
    filepoints, filecells = add_to_VtkMesh(filepoints, filecells, geometry)
    return filepoints, filecells
end

function to_VtkMesh(
    geometry :: Matrix{T}
    ) where T <: DiscreteGeometry3D
    return to_VtkMesh(vec(geometry))
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
function vtk_cell_type(a::BilinearQuad)
    return WriteVTK.VTKCellTypes.VTK_QUAD
end
function vtk_cell_type(a::BilinearQuadSurf)
    error("VTK has no equivalent to a BilinearQuadSurf. Please convert"*
        " to BilinearQuads. Eg: convert(Vector{UNSflow.BilinearQuad}"*
        ", ::BilinearQuadSurf)")
end
function vtk_cell_type(a::LatticeQuadSurf)
    error("VTK has no equivalent to a LatticeQuadSurf. Please convert"*
        " to Vector(Line2). Eg: convert(Vector{UNSflow.Line2}, "*
        "::LatticeQuadSurf). Alternately, your data may be better "*
        "represented as BilinearQuads.")
end

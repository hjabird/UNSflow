#===============================================================================
    UnstructuredMesh.jl

    Representation of an unstructured mesh with celldata & pointdata.

    Initial code: HJAB 2018

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


mutable struct UnstructuredMesh
    pointdata :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}}
    celldata  :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}}
    points :: Vector{Vector3D}
    cells  :: Vector{DiscreteGeometry3D}

    function UnstructuredMesh(
        pointdata :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}},
        celldata  :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}},
        points :: Vector{Vector3D},
        cells  :: Vector{DiscreteGeometry3D})
        new(pointdata, celldata, points, cells)
    end
end

function UnstructuredMesh()
    return UnstructuredMesh(
        Dict{String, Union{Vector{Float64}, Vector{Vector3D}}(),
        Dict{String, Union{Vector{Float64}, Vector{Vector3D}}(), 
        Vector{Vector3D}(), 
        Vector{DiscreteGeometry3D}())
end

"""
Add a cell to the mesh. Returns the cell's index and the indices of the
points accociated with it.
"""
function add_cell!(a::UnstructuredMesh, to_add::DiscreteGeometry3D)
    @assert(acceptable_cell_type(UnstructuredMesh, to_add),
        string("UNSflow.UnstructuredMesh sadly does not accept cells of type ",
        typeof(a), " to avoid problems with it having more than once peice",
        " of cell data. Try converting to Point3D, Line2, PolyLine2 or ",
        "BilinearQuad."))
    push!(a.cells, to_add)
    cell_idx = length(a.cells)
    coordinates = coords(a)
    num_coords = length(coordinates)
    append!(a.points, coordinates)
    point_idxs = collect(length(a.points) - num_coords + 1: length(a.points))
    return (cell_idx, point_idxs)
end

function add_celldata!(a::UnstructuredMesh, to_add_idx::Int, 
        data_name::String, value::Union{T, Vector3D}) where T <: Real
    
    # NB: MULTIPLE EXIT POINTS FROM FUNCTION
    @assert(0 < to_add_idx < length(a.cells), string("Index argument",
        " had value greater than the length of the cell vector. ",
        "0 < idx < length(cell vector). Argument was ", to_add_idx,
        " and length(::UnstructuredMesh.cells) was ", length(a.cells), "."))

    if haskey(a.celldata, data_name)
        if length(a.celldata[data_name]) > 0
            @assert(typeof(a.celldata[data_name][1]) == typeof(value),
                string("Tried to add data of type ", typeof(value),
                " to celldata with name ", data_name, " which contains data",
                " of type ", typeof(a.celldata[dataname][1])))
            if to_add_idx > length(a.celldata[data_name])
                length_deficit = length(a.cells) - length(a.celldata[data_name])
                append!(a.celldata[data_name], 
                    fill(typeof(value) == Vector3D ? 
                        Vector3D(NaN, NaN, NaN) : NaN, length_deficit))
            end
            a.celldata[data_name][to_add_idx] = value # Yay.    
            return # EXIT 1
        end
    end
    # Data of data_name doesn't exist or is zero length.
    if typeof(value) <: Real
        a.celldata[data_name] = fill(NaN, length(a.cells))
    else
        a.celldata[data_name] = fill(Vector3D(NaN, NaN, NaN), 
                                                            length(a.cells))
    end
    a.celldata[data_name][to_add_idx] = value
    return # EXIT 2
end

function add_pointdata!(a::UnstructuredMesh, point_idx::Int, 
    data_name::String, value::Union{Real, Vector3D})

    @assert(0 < point_idx < length(a.points), string("Index argument",
    " had value greater than the length of the points vector. ",
    "0 < idx < length(point vector). Argument was ", point_idx,
    " and length(::UnstructuredMesh.cells) was ", length(a.points), "."))

    if haskey(a.pointdata, data_name)
        if length(a.pointdata[data_name]) > 0
            @assert(typeof(a.pointdata[data_name][1]) == typeof(value),
                string("Tried to add data of type ", typeof(value),
                " to celldata with name ", data_name, " which contains data",
                " of type ", typeof(a.pointdata[dataname][1])))
            if point_idx > length(a.pointdata[data_name])
                length_deficit = length(a.cells) - length(a.pointdata[data_name])
                append!(a.pointdata[data_name], 
                    fill(typeof(value) == Vector3D ? 
                        Vector3D(NaN, NaN, NaN) : NaN, length_deficit))
            end
            a.pointdata[data_name][point_idx] = value # Yay.    
            return # EXIT 1
        end
    end
    # Data of data_name doesn't exist or is zero length.
    if typeof(value) <: Real
        a.pointdata[data_name] = fill(NaN, length(a.cells))
    else
        a.pointdata[data_name] = fill(Vector3D(NaN, NaN, NaN), 
                                                        length(a.cells))
    end
    a.pointdata[data_name][point_idx] = value
    return # EXIT 2
end


function acceptable_cell_type(Type{UnstructuredMesh}, a::DiscreteGeometry3D)
    typ = typeof(a)
    if typ == Point3D
        return true
    elseif typ == Line2
        return true
    elseif typ == PolyLine2
        return true
    elseif typ == BilinearQuad
        return true
    else
        return false
    end 
end

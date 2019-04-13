#===============================================================================
    UnstructuredMesh.jl

    Representation of an unstructured mesh with celldata & pointdata.
    It is expected that classes will add their own methods to add objects to 
    the mesh.

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
import WriteVTK

mutable struct UnstructuredMesh
    pointdata :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}}}
    celldata  :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}}}
    points :: Vector{Vector3D}
    cells  :: Vector{DiscreteGeometry3D}

    function UnstructuredMesh(
        pointdata :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}}},
        celldata  :: Dict{String, Union{Vector{Float64}, Vector{Vector3D}}},
        points :: Vector{Vector3D},
        cells  :: Vector{DiscreteGeometry3D})
        new(pointdata, celldata, points, cells)
    end
end

function UnstructuredMesh()
    return UnstructuredMesh(
        Dict{String, Union{Vector{Float64}, Vector{Vector3D}}}(),
        Dict{String, Union{Vector{Float64}, Vector{Vector3D}}}(), 
        Vector{Vector3D}(), 
        Vector{DiscreteGeometry3D}())
end

function UnstructuredMesh(cells::Any)
    ret = UnstructuredMesh()
    if typeof(cells) <: DiscreteGeometry3D
        add_cell!(ret, cells)
    elseif hasmethod(iterate, Tuple{typeof(cells)})
        for cell in cells
            add_cell!(ret, cell)
        end
    else
        error("Could not initialise Unstructured mesh with cells of type ", 
            typeof(cells), ". This should be a supertype of DiscreteGeometry3D",
            " or an iterable object.")
    end
    return ret
end

function Base.push!(a::UnstructuredMesh, 
    b::Any, 
    controldict=Dict{String, Any}())
    error("Type ", typeof(b), " has not provided a method to automatically ",
        "add it and its attributes to the UNSflow.UnstructuredMesh. You'll ",
        "have to use the add_cell!, add_pointdata!, add_celldata! based",
        " API.")
    return
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
    coordinates = coords(to_add)
    num_coords = length(coordinates)
    append!(a.points, coordinates)
    point_idxs = collect(length(a.points) - num_coords + 1: length(a.points))
    return (cell_idx, point_idxs)
end

function add_cells!(a::UnstructuredMesh, to_add::Any)
    ret = Vector{Tuple{Int64, Vector{Int64}}}()
    # Catch errors which'd probably just annoy the user if I pointed it out...
    if typeof(to_add) <: DiscreteGeometry3D
        push!(ret, add_cell!(a, to_add))
    else    
        @assert(hasmethod(iterate, Tuple{typeof(to_add)}),
            string("add_cell! takes arguement of ",
            "an iterable (ie. has method iterate) of DiscreteGeometry3D. ",
            "Called with add_cell!(::UnstructuredMesh, ::", typeof(to_add),
            ")."))
        @assert(hasmethod(length, Tuple{typeof(to_add)}), "The iterable argument ",
            "to_add must have a method length(::typeof(to_add)).")
        if length(to_add) > 0
            ret = map(x->add_cell!(a, x), to_add)
        end
    end
    return ret
end

function add_celldata!(a::UnstructuredMesh, to_add_idx::Int, 
        data_name::String, value::Union{T, Vector3D}) where T <: Real
    
    # NB: MULTIPLE EXIT POINTS FROM FUNCTION
    @assert(0 < to_add_idx <= length(a.cells), string("Index argument",
        " had value greater than the length of the cell vector. ",
        "0 < idx <= length(cell vector). Argument was ", to_add_idx,
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
    if !haskey(a.celldata, data_name) 
        assert_on_name_conflict(a, data_name)
    end
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

    @assert(0 < point_idx <= length(a.points), string("Index argument",
        " had value greater than the length of the points vector. ",
        "0 < idx <= length(point vector). Argument was ", point_idx,
        " and length(::UnstructuredMesh.cells) was ", length(a.points), "."))

    if haskey(a.pointdata, data_name)
        if length(a.pointdata[data_name]) > 0
            @assert(typeof(a.pointdata[data_name][1]) == typeof(value),
                string("Tried to add data of type ", typeof(value),
                " to celldata with name ", data_name, " which contains data",
                " of type ", typeof(a.pointdata[dataname][1])))
            if point_idx > length(a.pointdata[data_name])
                length_deficit = length(a.points)-length(a.pointdata[data_name])
                append!(a.pointdata[data_name], 
                    fill(typeof(value) == Vector3D ? 
                        Vector3D(NaN, NaN, NaN) : NaN, length_deficit))
            end
            a.pointdata[data_name][point_idx] = value # Yay.    
            return # EXIT 1
        end
    end
    # Data of data_name doesn't exist or is zero length.    
    if !haskey(a.pointdata, data_name) 
        assert_on_name_conflict(a, data_name)
    end
    if typeof(value) <: Real
        a.pointdata[data_name] = fill(NaN, length(a.cells))
    else
        a.pointdata[data_name] = fill(Vector3D(NaN, NaN, NaN), 
                                                        length(a.cells))
    end
    a.pointdata[data_name][point_idx] = value
    return # EXIT 2
end

function to_vtk_file(a::UnstructuredMesh, path::String)
    cells = Vector{WriteVTK.MeshCell}(undef, length(a.cells))
    # Because of the way we've directly stored the UNSflow.DiscreteGeometry3D
    # we have the pain of working out the indices of points...
    points = convert(Matrix{Float64}, a.points)
    for i = 1 : length(cells)
        pt_idxs = Vector{Int64}()
        coordinates = coords(a.cells[i])
        for j = 1: length(coordinates)
            c = findfirst(x->x === coordinates[j], a.points)
            push!(pt_idxs, c)
        end
        cells[i] = WriteVTK.MeshCell(vtk_cell_type(a.cells[i]), pt_idxs)
    end
    vtkfile = WriteVTK.vtk_grid(path, points, cells)
    # Add celldata and pointdata.
    for dataset in a.celldata
        dataname = dataset[1]
        data = dataset[2]
        if length(data) < length(a.cells)
            append!(data, fill(typeof(data[1]) == Vector3D ? 
                Vector3D(NaN, NaN, NaN) : Float64(NaN),
                length(a.cells) - length(data)))
        end
        if typeof(data) == Vector{Vector3D}
            data = convert(Matrix{Float64}, data)
        end
        vtkfile = WriteVTK.vtk_cell_data(vtkfile, data, dataname)
    end
    for dataset in a.pointdata
        dataname = dataset[1]
        data = dataset[2]
        if length(data) < length(a.points)
            append!(data, fill(typeof(data[1]) == Vector3D ? 
                Vector3D(NaN, NaN, NaN) : Float64(NaN),
                length(a.cells) - length(data)))
        end
        if typeof(data) == Vector{Vector3D}
            data = convert(Matrix{Float64}, data)
        end
        vtkfile = WriteVTK.vtk_cell_data(vtkfile, data, dataname)
    end
    outfile = WriteVTK.vtk_save(vtkfile)
    return outfile
end

function acceptable_cell_type(::Type{UnstructuredMesh}, a::DiscreteGeometry3D)
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

function grow_field_vectors!(a::UnstructuredMesh)
    lenc = length(a.cells)
    lenp = length(a.points)

    for field in a.pointdata
        deficit = lenp - length(field[2])
        if deficit > 0
            append!(field[2], fill(typeof(field[2][1]) == Vector3D ? 
                Vector3D(NaN, NaN, NaN) : NaN, deficit))
        end
    end
    for field in a.celldata
        deficit = lenc - length(field[2])
        if deficit > 0
            append!(field[2], fill(typeof(field[2][1]) == Vector3D ? 
                Vector3D(NaN, NaN, NaN) : NaN, deficit))
        end
    end
end

function assert_on_name_conflict(a::UnstructuredMesh, new_name::String)
    nc = check_name_conflict(a, new_name)
    @assert( nc == false ,
        string("Name ", new_name, " was already in use in either pointdata or ",
            "celldata. You cannot use the same name in both - it'll crash ",
            "Paraview (home of the most fragile file readers around)."))
    return
end

function check_name_conflict(a::UnstructuredMesh, new_name::String)
    if new_name in keys(a.pointdata)
        return true
    elseif new_name in keys(a.celldata)
        return true
    else
        return false
    end
end

function merge_redundant_points!(a::UnstructuredMesh)
    grow_field_vectors!(a::UnstructuredMesh)
    # Some points might be repeated, but if they have different pointdata 
    # for two points at the same location, we don't want to merge that.
    # We'll try and keep the one with the first index...
    # Step 1 find repeated points and record index of first instance.
    firstidx = Vector{Int64}(undef, length(a.points))
    encountered = Dict{Vector3D, Int64}()
    for i = 1 : length(a.points)
        if haskey(encountered, a.points[i])
            firstidx[i] = encountered[a.points[i]]
        else
            firstidx[i] = i
            encountered[a.points[i]] = i
        end
    end
    # Step 2 correct where points have different data associated with them.
    for field in a.pointdata
        for i = 1 : length(a.points)
            firstidx[i] = field[2][i] == field[2][firstidx[i]] ? firstidx[i] : i
        end
    end
    # Step 3 adjust firstidx for where the offset is changing...
    offset = 0
    for i = 1 : length(a.points)
        firstidx[i] = firstidx[i] - offset
        if firstidx[i] != i
            offset += 1
        end
    end
    # Step 4 adjust cell point indexes & make new point vector
    new_points = Vector{Int64}()
    offset = 1
    for i = 1 : length(a.points)
        if firstidx[i] >= offset
            offset += 1
            push!(new_points, firstidx)
        end
    end
    for 
    error("CONSTRUCTION ZONE: KEEP OUT")
    return
end

#===============================================================================
    MeshDataLinker.jl

    Allow data to be stored externally from the objects used for computation
    whilst also making it easy to add to an unstructured mesh.
    
    add_pointdata!(::MeshDataLinker, ::obj, dataname, data) and 
    add_celldata!(::MeshDataLinker, ::obj, dataname, data) need to be redefined 
    on a per-object basis for easy use.

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

mutable struct MeshDataLinker
    pointdata ::Dict{String, Union{ Dict{Vector3D, Float64 },
                                    Dict{Vector3D, Vector3D}}}
    celldata :: Dict{String, Union{ Dict{DiscreteGeometry3D, Float64 },
                                    Dict{DiscreteGeometry3D, Vector3D}}}
end

function MeshDataLinker()
    return MeshDataLinker(
        Dict{   String, 
                Union{  Dict{Vector3D, Float64 }, 
                        Dict{Vector3D, Vector3D}}}(),
        Dict{   String, 
                Union{  Dict{DiscreteGeometry3D, Float64 },
                        Dict{DiscreteGeometry3D, Vector3D}}}())
end

function add_celldata!(a::MeshDataLinker, dataname::String, 
    cellgeom::DiscreteGeometry3D, value::Union{Float64, Vector3D}; 
    warn_value_overwrite::Bool=false)

    if warn_value_overwrite && haskey(a.celldata, dataname) &&
        haskey(a.celldata[dataname], cellgeom)
        @warn "Overwriting celldata."
    end
    if !haskey(a.celldata, dataname)
        a.celldata[dataname] = Dict{DiscreteGeometry3D, typeof(value)}()
    end
    a.celldata[dataname][cellgeom] = value
end

function add_pointdata!(a::MeshDataLinker, dataname::String, 
    point::Vector3D, value::Union{Float64, Vector3D}; 
    warn_value_overwrite::Bool=false)
    
    if warn_value_overwrite && haskey(a.celldata, dataname) &&
        haskey(a.pointdata[dataname], point)
        @warn "Overwriting pointdata."
    end
    if !haskey(a.celldata, dataname)
        a.pointdata[dataname] = Dict{Vector3D, typeof(value)}()
    end
    a.pointdata[dataname][point] = value
end

function add_data!(a::UnstructuredMesh, b::MeshDataLinker)
    grow_field_vectors!(a)
    # Points:
    point_idxs = Dict{Vector3D, Int64}()
    for i = 1 : length(a.points)
        point_idxs[a.points[i]] = i
    end
    for field in b.pointdata
        for pnt in field[2]
            if haskey(point_idxs, pnt[1])
                a.pointdata[field[1]][point_idxs[pnt[1]]] = pnt[2]
            end
        end
    end
    # Geometry
    geometry_idxs = Dict{DiscreteGeometry3D, Int64}()
    for i = 1 : length(a.cells)
        geometry_idxs[a.cells[i]] = i
    end
    for field in b.celldata
        for cell in field[2]
            if haskey(geometry_idxs, cell[1])
                a.celldata[field[1]][geometry_idxs[cell[1]]] = cell[2]
            end
        end
    end
    return
end

#===============================================================================
    Vorticity3DCollector.jl

    Represents a collection of Vorticity3D(ies) that may themselves be
    Vorticity3DCollectors.

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
include("Vorticity3D.jl")
include("Vector3D.jl")

abstract type Vorticity3DCollector <: Vorticity3D
end

#=
The interface:

Vorticity3DCollector.children
    are an iterable of the children of the collector.

child_vorticities(a::Vorticity3DCollector)
    returns a vector of the vorticity of all the non-Vorticity3DCollector
    Vorticity3D objects owned by this collector, and, recursively, its
    children.
=#

function get_children(
    a::Vorticity3DCollector,
    typefilter=Vorticity3D)

    retv = Vector{typefilter}()
    for child in a
        if typeof(child) <: typefilter
            push!(retv, child)
        end
    end
    return retv
end

function get_children_recursive(
    a::Vorticity3DCollector,
    typefilter=Vorticity3D
    )
    retv = Vector{typefilter}()
    for child in a
        if typeof(child) <: Vorticity3DCollector
            rv = get_children_recursive(child, typefilter)
            retv = vcat(retv, rv)
        elseif typeof(child) <: typefilter
            push!(retv, child)
        end
    end
    return retv
end

#= Defualt container interaction methods -------------------------------------=#
function Base.push!(a::Vorticity3DCollector, to_be_added::Vorticity3D)
    push!(a.children, to_be_added)
end

function Base.append!(
    a::Vorticity3DCollector,
    iterable_of_things_to_be_appended)

    @assert(eltype(iterable_of_things_to_be_appended) <: Vorticity3D,
        "To add individual Vorticity3D objects you want to use " *
        "push!(a::Vorticity3DCollector, b::Vorticity3D). To add the children " *
        "of a Vorticity3DCollector (\"b\") to this Vorticity3DCollector try " *
        "append!(a::Vorticity3DCollector, get_children(b)).")
    append!(a.children, iterable_of_things_to_be_appended)
end

function Base.pop!(a::Vorticity3DCollector)
    return pop!(a.children)
end

function Base.isempty(a::Vorticity3DCollector)
    return isempty(a.children)
end

function Base.empty!(a::Vorticity3DCollector)
    empty!(a.children)
end


#= Default vorticity body methods --------------------------------------------=#
function vorticity(a::Vorticity3DCollector)
    return mapreduce(vorticity, +, a.children, init=Vector3D(0,0,0))
end

function centre(a::Vorticity3DCollector)
    vort = 0.0
    center = Vector3D(0,0,0)    # American spelling has its uses.
    for child in a.children
        vort += vorticity(child)
        center += centre(child) * vorticity
    end
    center /= vorticity
    return center
end

function effective_radius(a::Vorticity3DCollector)
    n_children = length(a.children)
    radii = zeros(Float64, n_children)
    c = centre(a)
    for i = 1:n_children
        radii[i] = abs(centre(a.children[i]) - c) +
            effective_radius(a.children[i])
    end
    rad = maximum(radii)
    return rad
end

function euler!(
    a::Vorticity3DCollector,
    b::Vorticity3D,
    dt::Real)
    # So we can call the update methods of our children...
    cpy = deepcopy(a)
    for child in cpy.children
        euler!(child, b, dt)
    end
    a.children = cpy.children
    return
end

function induced_velocity(
    a::Vorticity3DCollector,
    measurement_point :: Vector3D
    )
    vel = Vector3D(0,0,0)
    for child in a
        vel += induced_velocity(child, measurement_point)
    end
    return vel
end

function induced_velocity_curl(
    a::Vorticity3DCollector,
    measurement_point :: Vector3D
    )
    curly = zeros(3,3)
    for child in a
        curly += induced_velocity_curl(child, measurement_point)
    end
    return curly
end

#= Default iterator ----------------------------------------------------------=#
function Base.getindex(
    a::Vorticity3DCollector,
    i::Integer)
    return a.children[i]
end

function Base.setindex!(
    a::Vorticity3DCollector,
    i::Integer,
    v::Vorticity3D)
    a.children[i] = v
    return
end

function Base.length(a::Vorticity3DCollector)
    return length(a.children)
end

function Base.size(a::Vorticity3DCollector)
    return size(a.children)
end

if VERSION >= VersionNumber(0, 7, 0)    # Why not use Julia you ask?
    function Base.firstindex(a::Vorticity3DCollector)
        return 1
    end

    function Base.lastindex(a::Vorticity3DCollector)
        return length(a.children)
    end

    function Base.iterate(a::Vorticity3DCollector, state=1)
        if state > length(a)
            return nothing
        else
            return (a.children[state], state + 1)
        end
    end
else # Version 0.6 or older.
    function Base.start(a::Vorticity3DCollector)
        return 1
    end

    function Base.next(a::Vorticity3DCollector, state::Integer)
        return (a[state], state+1)
    end

    function Base.done(a::Vorticity3DCollector, state::Integer)
        return state > length(a)
    end
end

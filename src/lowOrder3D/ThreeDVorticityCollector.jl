#===============================================================================
    ThreeDVorticityCollector.jl

    Represents a collection of ThreeDVorticity(ies) that may themselves be
    ThreeDVorticityCollectors.

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
include("ThreeDVorticity.jl")
include("ThreeDVector.jl")

abstract type ThreeDVorticityCollector <: ThreeDVorticity
end

#=
The interface:

ThreeDVorticityCollector.children
    are an iterable of the children of the collector.

child_vorticities(a::ThreeDVorticityCollector)
    returns a vector of the vorticity of all the non-ThreeDVorticityCollector
    ThreeDVorticity objects owned by this collector, and, recursively, its
    children.
=#

function get_children(
    a::ThreeDVorticityCollector,
    typefilter=ThreeDVorticity)

    retv = Vector{typefilter}()
    for child in a
        if typeof(child) <: typefilter
            push!(retv, child)
        end
    end
    return retv
end

function get_children_recursive(
    a::ThreeDVorticityCollector,
    typefilter=ThreeDVorticity
    )
    retv = Vector{typefilter}()
    for child in a
        if typeof(child) <: ThreeDVorticityCollector
            rv = get_children_recursive(child, typefilter)
            vcat(revt, rv)
        elseif typeof(child) <: typefilter
            push!(retv, child)
        end
    end
    return retv
end

#= Default vorticity body methods --------------------------------------------=#
function vorticity(a::ThreeDVorticityCollector)
    return mapreduce(vorticity, 0.0, +, a.children)
end

function centre(a::ThreeDVorticityCollector)
    vort = 0.0
    center = ThreeDVector(0,0,0)    # American spelling has its uses.
    for child in a.children
        vort += vorticity(child)
        center += centre(child) * vorticity
    end
    center /= vorticity
    return center
end

function effective_radius(a::ThreeDVorticityCollector)
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
    a::ThreeDVorticityCollector,
    b::ThreeDVorticity,
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
    a::ThreeDVorticityCollector,
    measurement_point :: ThreeDVector
    )
    vel = @parallel (+) for child in a
        vel += induced_velocity(child, measurement_point)
    end
    return vel
end

function induced_velocity_curl(
    a::ThreeDVorticityCollector,
    measurement_point :: ThreeDVector
    )
    curly = @parallel (+) for child in a
        curly += induced_velocity_curl(child, measurement_point)
    end
    return curly
end

#= Default iterator ----------------------------------------------------------=#
function Base.getindex(
    a::ThreeDVorticityCollector,
    i::Integer)
    return a.children[i]
end

function Base.setindex!(
    a::ThreeDVorticityCollector,
    i::Integer,
    v::ThreeDVorticity)
    a.children[i] = v
    return
end

function Base.length(a::ThreeDVorticityCollector)
    return length(a.children)
end

function Base.size(a::ThreeDVorticityCollector)
    return size(a.children)
end

if VERSION >= VersionNumber(0, 7, 0)    # Why not use Julia you ask?
    function Base.firstindex(a::ThreeDVorticityCollector)
        return 1
    end

    function Base.lastindex(a::ThreeDVorticityCollector)
        return length(a.children)
    end

    function Base.iterate(a::ThreeDVorticityCollector, state=1)
        if state > length(a)
            return nothing
        else
            return (a.children[state], state + 1)
        end
    end
else # Version 0.6 or older.
    function Base.start(a::ThreeDVorticityCollector)
        return 1
    end

    function Base.next(a::ThreeDVorticityCollector, state::Integer)
        return (a[state], state+1)
    end

    function Base.done(a::ThreeDVorticityCollector, state::Integer)
        return state > length(a)
    end
end

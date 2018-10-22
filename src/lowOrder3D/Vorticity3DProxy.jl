#===============================================================================
    Vorticity3DProxy.jl

    A proxy for a Vorticity3D object that, when subtyped, allows the 
    interception of function calls whilst providing a transparent interface
    to the source object's fields (so long as they don't have the same names
    as the proxy's field).

    Examples of use include LinkedVorticity3DProxy.

    WHAT THIS ACTUALLY DOES:
    It automatically directs reading and assignment of fields to the source 
    object. So:
    julia> mutable struct a <: Vorticity3D
    julia>     x :: Int64
    julia> end
    julia> mutable struct b <: Vorticity3DProxy
    julia>     source_object :: a 
    julia> end
    julia> c = a(4) # Make our really simple "Vorticity3D" object
    julia> c.x
    4
    julia> c.x = 6
    6
    julia> d = b(c) # Make a proxy to the Vorticity3D object
    julia> d.x      # We automatically get c's field names.
    6
    julia> d.x = 10 # We can even rebind them
    10
    julia> c.x # And the underlying data gets rebound too!
    10

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

abstract type Vorticity3DProxy <: Vorticity3D
end

# The magic bit ----------------------------------------------------------------
function Base.getproperty(a::Vorticity3DProxy, name::Symbol)
    if name in fieldnames(typeof(a))
        getfield(a, name)
    else
        getproperty(a.source_object, name)
    end
end

function Base.setproperty!(a::Vorticity3DProxy, name::Symbol, value)
    if name in fieldnames(typeof(a))
        setfield!(a, name)
    else
        setproperty!(a.source_object, name, value)
    end
end

# The mundane automatic redirection of Vorticity3D calls... --------------------
function centre(a::Vorticity3DProxy)
    return centre(a.source_object)
end

function effective_radius(a::Vorticity3DProxy)
    return effective_radius(a.source_object)
end

function vorticity(a::Vorticity3DProxy)
    return vorticity(a.source_object)
end

function induced_velocity(a::Vorticity3DProxy, measurement_point::Vector3D)
    return induced_velocity(a.source_object, measurement_point)
end

function induced_velocity_curl(
    a::Vorticity3DProxy, 
    measurement_point::Vector3D)
    return induced_velocity_curl(a.source_object, measurement_point)
end

function euler!(a::Vorticity3DProxy, influence_field::Vorticity3D, dt::Real)
    return euler!(a.source_object, influence_field, dt)
end

function state_vector_length(a::Vorticity3DProxy)
    return state_vector_length(a.source_object)
end

function state_vector(a::Vorticity3DProxy)
    return state_vector(a.source_object)
end

function update_using_state_vector!(
    a::Vorticity3DProxy,
    state_vect::Vector{T}) where T <: Real
    return update_using_state_vector!(a.source_object, state_vect)
end

function state_time_derivative(
    a::Vorticity3DProxy,
    inducing_bodies::Vorticity3D)

    return state_time_derivative(a.source_object, inducing_bodies)
end

function vorticity_vector_length(a::Vorticity3DProxy)
    return vorticity_vector_length(a.source_object)
end

function vorticity_vector(a::Vorticity3DProxy)
    return vorticity_vector(a.source_object)
end

function update_using_vorticity_vector!(
    a::Vorticity3DProxy,
    vort_vect::Vector{T}) where T <: Real
    
    update_using_vorticity_vector!(a.source_object, vort_vect)
    return
end

function vorticity_vector_velocity_influence(
    a::Vorticity3DProxy,
    mes_pnt::Vector3D
    )

    return vorticity_vector_velocity_influence(a.source_object, mes_pnt)
end

function steady_forces(a::Vorticity3DProxy, args...; kwarg...)
    args = tuple(a.source_object, args...)
    argt = tuple(map(typeof, args)...)
    if hasmethod(steady_forces, argt)
        ret = steady_forces(args...; kwarg...)
    else
        error(string(typeof(a), " has no method steady_forces(...) with",
            " arguments of type (contained in tuple here):", argt))
    end
    return ret
end

function steady_pressures(a::Vorticity3DProxy, args...; kwarg...)
    args = tuple(a.source_object, args...)
    argt = tuple(map(typeof, args)...)
    if hasmethod(steady_pressures, argt)
        ret = steady_pressures(args...; kwarg...)
    else
        error(string(typeof(a), " has no method steady_pressures(...) with",
            " arguments of type (contained in tuple here):", argt))
    end
    return ret
end

function steady_loads(a::Vorticity3DProxy, args...; kwarg...)
    args = tuple(a.source_object, args...)
    argt = tuple(map(typeof, args)...)
    if hasmethod(steady_loads, argt)
        ret = steady_loads(args...; kwarg...)
    else
        error(string(typeof(a), " has no method steady_loads(...) with",
            " arguments of type (contained in tuple here):", argt))
    end
    return ret
end

# Iterator / Vorticity3DCollector functions ------------------------------------
# Duck typing: adding in appropriate functions just incase my source is of a 
# type, when actually, it might not be. If it has legs like a duck, wings like
# a duck and a bill like a duck, but can't waddle fly or quack, is it a duck!?
function Base.getindex(
    a::Vorticity3DProxy,
    i::Integer)    
    if hasmethod(Base.setindex!, 
        Tuple{typeof(a.source_object), Integer})
        return Base.getindex(a, i)    
    else
        error(string("Source object of type ", typeof(a.source_object),
            "has no method getindex."))
        return
    end
end

function Base.setindex!(
    a::Vorticity3DProxy,
    i::Integer,
    v::Vorticity3D)   
    if hasmethod(Base.setindex!, 
        Tuple{typeof(a.source_object), Integer, Vorticity3D})
        return Base.setindex!(a, i, v)    
    else
        error(string("Source object of type ", typeof(a.source_object),
            "has no method setindex!."))
    end
    return
end

function Base.length(a::Vorticity3DProxy)
    if hasmethod(Base.length, Tuple{typeof(a.source_object)})
        return length(a.source_object)    
    else
        error(string("Source object of type ", typeof(a.source_object),
            "has no method length."))
    end
end

function Base.size(a::Vorticity3DProxy)
    if hasmethod(Base.size, Tuple{typeof(a.source_object)})
        return size(a.source_object)    
    else
        error(string("Source object of type ", typeof(a.source_object),
            "has no method size."))
    end
end

function Base.firstindex(a::Vorticity3DProxy)
    if hasmethod(Base.firstindex, Tuple{typeof(a.source_object)})
        return Base.firstindex(a.source_object)
    else
        error(string("Source object of type ", typeof(a.source_object),
            "has no method firstindex. Are you sure it is iterable?"))
    end
end

function Base.lastindex(a::Vorticity3DProxy)
    if hasmethod(Base.lastined, Tuple{typeof(a.source_object)})
    return lastindex(a.source_object)
else
    error(string("Source object of type ", typeof(a.source_object),
        "has no method lastindex. Are you sure it is iterable?"))
end
end

function Base.iterate(a::Vorticity3DProxy, state=1)
    if hasmethod(Base.firstindex, Tuple{typeof(a.source_object), Any})
        iterate(a.source_object, state=state)
    else
        error(string("Source object of type ", typeof(a.source_object),
            "has no method iterate. Are you sure it is iterable?"))
    end
end

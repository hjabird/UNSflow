#===============================================================================
    Vorticity3D.jl

    Represents some kind of vorticity within the flow domain. Does not define
    goemetry or any other such thing. May be subtyped by groups of vortex
    objects or by individual objects.

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

abstract type Vorticity3D
end

#= Expected interface for Vorticity3D:

Constructor is defined by the concrete type.
=#

"""
returns something somehow representative of the bodies' centre location
as a Vector3D. Whether this is geometrically central, or vorticity
weighted does not matter so long as it is considered by the programmer along
with the effective_radius definition.
"""
function centre(this::Vorticity3D)
    error("Not implemented for type ", typeof(this), ".")
    return Vector3D(0,0,0)
end

"""
returns a Real representative of the effective radius from the centre
of a sphere that contains the 3D vorticity body's vorticy.
"""
function effective_radius(this::Vorticity3D)
    error("Not implemented for type ", typeof(this), ".")
    return Float64(0.0)
end

"returns the integral of vorticity within the body as Vector3D."
function vorticity(this::Vorticity3D)
    error("Not implemented for type ", typeof(this), ".")
    return Vector3D(0,0,0)
end

"""
returns the velocity induced by the body at measurement_point as a
Vector3D.
"""
function induced_velocity(this::Vorticity3D, measurement_point::Vector3D)
    error("Not implemented for type ", typeof(this), ".")
    return Vector3D(0,0,0)
end

"""
returns the curl in the velocity induced by the body at measurement_point as
a 3 by 3 matrix, by which vorticity at that point can be multiplied.
"""
function induced_velocity_curl(this::Vorticity3D, measurement_point::Vector3D)
    error("Not implemented for type ", typeof(this), ".")
    return zeros(3,3)
end

"""
Applies convection to this body due to velocities induced by
influence_field. Uses the forward Euler method.
"""
function euler!(this::Vorticity3D, influence_field::Vorticity3D, dt::Real)
    state = state_vector(this)
    dstate = state_time_derivative(this, inducing_field)
    state += dstate * dt
    update_using_state_vector(this, state)
    return
end

"""
Returns the length of the state vector associated with a vorticity object.
For solving convection using an ODE solver
"""
function state_vector_length(a::Vorticity3D)
    error("state_vector_length(a::Vorticity3D) was not implemented for ",
        "subtype ", typeof(this), ".")
    return
end

"""
Returns the "state vector" representing the state of the Vorticity3D object
that might change due to convection.

This is intended for use during the solution of the ODE problem of convection.

Associated useful functions are state_vector_length(::Vorticity3D),
update_using_state_vector!(::Vorticity3D, state_vector::Vector{Float64})
and state_vector_time_derivative(this::Vorticity3D, 
inducing_bodies::Vorticity3D).
"""
function state_vector(a::Vorticity3D)
    @assert(a <: Vorticity3DCollector, "The specialised function should have "*
        "been called here. Check Vorticity3DCollector.jl") 
    state_vect = Vector{Float64}(undef, state_vector_length(a))
    coord_vect = coords(a.geometry)
    for c in coord_vect
        append!(state_vect, convert(Vector{Float64}, c))
    end
    return state_vect
end

"""
Update the Vorticity3D object using a state vector, usually generated
using some form of ODE solver.
"""
function update_using_state_vector!(
    this::Vorticity3D,
    state_vect::Vector{Float64})
    @assert(typeof(this) != VortexParticle3D, "Wrong function dispatched!?")
    @assert(!(typeof(this) <: Vorticity3DCollector), 
        "Wrong function dispatched!?")
    
    coord_vect = coords(this.geometry)
    @assert(length(state_vect) == 3 * length(coord_vect), string(
        "Input state vector was the incorrect length. Length was ",
        length(state_vect), " but should have been ",
        3 * length(coord_vect), "."))
    for i = 1 : length(coord_vect)
        coord_vect[i].x = state_vect[i * 3 - 2]
        coord_vect[i].y = state_vect[i * 3 - 1]
        coord_vect[i].z = state_vect[i * 3]
    end
    coord_vect = coord_vect[length(coord_vect) * 3 + 1: end]
    return
end

"""
Returns the time derivative of the state vector of `this` due to
an inducing Vorticity3D.
"""
function state_time_derivative(
    this::Vorticity3D,
    inducing_bodies::Vorticity3D)

    @assert(typeof(this) != VortexParticle3D)
    @assert(!(typeof(this) <: Vorticity3DCollector))
    deriv_vect = Vector{Float64}()
    coord_vect = coords(this.geometry)
    for c in coord_vect
        append!(deriv_vect,
            convert(Vector{Float64}, induced_velocity(this, inducing_bodies)))
    end
    return deriv_vect
end

"""
Returns the length of a vector representing the vorticity controls within
the Vorticity3D
"""
function vorticity_vector_length(this::Vorticity3D)
    error("vorticity_vector_length(a::Vorticity3D) was not implemented for ",
        "subtype ", typeof(this), ".")
    return
end

"""
Returns a vector that represents the vorticity state of a vorticity3D
"""
function vorticity_vector(this::Vorticity3D)
    error("vorticity_vector(a::Vorticity3D) was not implemented for ",
        "subtype ", typeof(this), ".")
    return
end

"""
Update the vorticities within a Vorticity3D using a vector.
"""
function update_using_vorticity_vector!(
    this::Vorticity3D,
    vort_vect::Vector{Float64})

    error("update_using_vorticity_vector!(a::Vorticity3D, ",
        "vort_vect::Vector{Float64}) has not been reimlemented for ",
        "Vorticity3D subtype ", typeof(this),
        ".")
    return
end

"""
Compute the influence matrix of a vorticity3D on the velocity at a point in
space.

For a Vorticity3D with n vorticity controls (can be found using
vorticity_vector_length(this::Vorticity3D)), a 3 x n matrix is returned.
Multiplying this matrix by the vorticity vector will result in a vector
of length 3 representing the velocity at a point.
"""
function vorticity_vector_velocity_influence(
    this::Vorticity3D,
    mes_pnt::Vector3D
    )
    error("vorticity_vector_velocity_influence(a::Vorticity3D, ",
        "mes_pnt::Vector{Float64}) has not been reimlemented for ",
        "Vorticity3D subtype ", typeof(this),
        ".")
    return
end

"""
Generate an matrix describing how the vorticity_vector changes the normal 
velocity at given points and normals.
"""
function normal_velocity_influence_matrix(
    inducing_body::Vorticity3D, 
    mes_points::Vector{Vector3D},
    mes_normals::Vector{Vector3D})
    @assert(size(mes_normals) == size(mes_points))

    @assert(all(isfinite.(mes_normals)), "All input normals should have"*
        " only finite components.")
    if !all(isfinite.(mes_points))
        @warn "Points with NaN or Inf components."
    end

    # Call number of points n
    # Call vorticity_vector_length(inducing_body) m 

    # Expect n long vector of vector{Vector3D}(m)
    ivs = map(
        x->vorticity_vector_velocity_influence(inducing_body, x), 
        mes_points)
    # Expect n long vector of vector{Float64}(m)
    ins = map(x->dot(Matrix{Float64}(x[1]'), x[2]), zip(ivs, mes_normals))
    a = zeros(length(mes_points), vorticity_vector_length(inducing_body))
    for i = 1 : length(mes_points)
        a[i, :] = ins[i]
    end
    return a
end

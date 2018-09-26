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

centre(this::Vorticity3D)
    returns something somehow representative of the bodies' centre location
    as a Vector3D. Whether this is geometrically central, or vorticity
    weighted does not matter so long as it is considered by the programmer along
    with the effective_radius definition.

effective_radius(this::Vorticity3D)
    returns a Real representative of the effective radius from the centre
    of a sphere that contains the 3D vorticity body's vorticy.

vorticity(this::Vorticity3D)
    returns the integral of vorticity within the body as Vector3D.

induced_velocity(this::Vorticity3D, measurement_point::Vector3D)
    returns the velocity induced by the body at measurement_point as a
    Vector3D.

induced_velocity_curl(this::Vorticity3D, measurement_point::Vector3D)
    returns the curl in the velocity induced by the body at measurement_point as
    a 3 by 3 matrix, by which vorticity at that point can be multiplied.
=#

# Applies convection to this body due to velocities induced by
# influence_field. Uses the forward Euler method.
function euler!(this::Vorticity3D, influence_field::Vorticity3D, dt::Real)
    coords = coords(this.goemetry)
    vel = Vector{Vector3D}(undef, length(coords))
    for pair in zip(vel, coords)
        pair[1] = map(x->induced_velocity(influence_field, x), pair[2])
    end
    coords += vel
    return
end

function state_vector(a::Vorticity3D)
    # We assume here that the change in due only due to convection of points.
    # Otherwise we need to specially define this such as for VortexParticle3D
    state_vect = Vector{Float64}()
    if typeof(a) <: Vorticity3DCollector
        for child in get_children(a)
            append!(state_vect, state_vector(child))
        end
    else
        coord_vect = coords(a.geometry)
        for c in coord_vect
            append!(state_vect, convert(Vector{Float64}, c))
        end
    end
    return state_vect
end

function update_using_state_vector!(
    this::Vorticity3D,
    state_vect::Vector{Float64})
    # We assume here that the change in due only due to convection of points.
    # Otherwise we need to specially define this such as for VortexParticle3D
    if typeof(this) <: Vorticity3DCollector
        for child in get_children(this)
            state_vect = update_using_state_vector(child, state_vect)
        end
    else
        coord_vect = coords(this.geometry)
        for i = 1 : length(coord_vect)
            coord_vect[i].x = state_vect[i * 3 - 2]
            coord_vect[i].y = state_vect[i * 3 - 1]
            coord_vect[i].z = state_vect[i * 3]
        end
        coord_vect = coord_vect[length(coord_vect) * 3 + 1: end]
    end
    return state_vect
end

function state_time_derivative(
    this::Vorticity3D,
    inducing_bodies::Vorticity3D)
    # We assume here that the change in due only due to convection of points.
    # Otherwise we need to specially define this such as for VortexParticle3D
    deriv_vect = Vector{Float64}()
    if typeof(this) <: Vorticity3DCollector
        for child in get_children(this)
            append!(deriv_vect, state_time_derivative(child, inducing_bodies))
        end
    else
        coord_vect = coords(this.geometry)
        for c in coord_vect
            append!(deriv_vect,
                convert(Vector{Float64}, induced_velocity(this, inducing_bodies)))
        end
    end
    return deriv_vect
end

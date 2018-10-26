#===============================================================================
    FreeStream3D.jl

    A free stream representation that can be used within the Vorticity3D
    framework

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

mutable struct FreeStream3D <: Vorticity3D
    velocity :: Vector3D
end

function vorticity(a::FreeStream)
    return Vector3D(0,0,0)
end

function induced_velocity(a::FreeStream)
    return velocity
end

function induced_velocity_curl(a::FreeStream)
    return zeros(3,3)
end

function euler!(a::FreeStream, influence_field::Vorticity3D, dt::Real)
    return
end

function state_vector_length(a::FreeStream)
    return 0
end

function update_using_state_vector!(
    this::FreeStream,
    state_vect::Vector{T}) where T <: Real
    @assert(length(state_vect) == 0, "FreeStream has no state!")
    return
end

function state_time_derivative(
    a::FreeStream,
    inducing_bodies::Vorticity3D)
    return Vector{Float64}()
end

function vorticity_vector_length(a::FreeStream)
    # In some senses this is cheating, but this is needed in the
    # sense that one ought to be able to compute velocity from a
    # strength / influence multiplication.
    return 3
end

function vorticity_vector(a::FreeStream)
    return Vector{Float64}(a.velocity)
end

function update_using_vorticity_vector!(
    this::FreeStream,
    vort_vect::Vector{T}) where T <: Real
    @assert(length(vort_vect) == 3)
    a.velocity == convert(Vector3D, vort_vect)
    return
end

function vorticity_vector_velocity_influence(
    this::Vorticity3D,
    mes_pnt::Vector3D
    )
    ret = LinearAlgebra.Identity * a.velocity
    return ret
end

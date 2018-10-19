#===============================================================================
    LinkedVorticity3D.jl

    Links the vorticity of a child object to that of a parent object. All
    call refer only to the influence of the child object. Ie: One must 
    add vorticity_vector_velocity_influence(this) + 
    vorticity_vector_velocity_influence(parent).

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

mutable struct LinkedVorticity3D <: Vorticity3D
    child_object :: Vorticity3D
    source_object :: Vorticity3D
    vorticity_vector_transform_matrix :: Matrix{Float64}
    geometry # Links to child's geometry

    function LinkedVorticity3D(
        child_object::Vorticity3D, 
        source_object::Vorticity3D, 
        transformation_matrix::Matrix{T}) where T <: Real

        tmat = transformation_matrix
        @assert(size(tmat, 1) == vorticity_vector_length(child_object),
            string("LinkedVorticity3D: the vorticity vector length of the ",
            "child object does not match the number or rows in the ",
            "transformation matrix. Child object vorticity vector length is ",
            vorticity_vector_length(child_object), " and size of the  ",
            "matrix is ", size(tmat)))
        @assert(size(tmat, 2) == vorticity_vector_length(source_object),
            string("LinkedVorticity3D: the vorticity vector length of the ",
            "source object does not match the number or columns in the ",
            "transformation matrix. Source object vorticity vector length is ",
            vorticity_vector_length(source_object), " and size of the  ",
            "matrix is ", size(tmat)))
        if any(isnan.(tmat)) || any(isinf.(tmat))
            @warn "Transformation matrix contained NaN or Inf values."
        end
        geom = child_object.geometry
        new(child_object, source_object, transformation_matrix, geom)
    end
end

function centre(a::LinkedVorticity3D)
    return centre(a.child_object)
end

function effective_radius(a::LinkedVorticity3D)
    return effective_radius(a.child_object)
end

function vorticity(a::LinkedVorticity3D)
    return vorticity(a.child_object)
end

function induced_velocity(a::LinkedVorticity3D, measurement_point::Vector3D)
    return induced_velocity(a.child_object, measurement_point)
end

function induced_velocity_curl(
    a::LinkedVorticity3D, 
    measurement_point::Vector3D)
    return induced_velocity_curl(a.child_object, measurement_point)
end

function euler!(a::LinkedVorticity3D, influence_field::Vorticity3D, dt::Real)
    return euler!(a.child_object, influence_field, dt)
end

function state_vector_length(a::LinkedVorticity3D)
    return state_vector_length(a.child_object)
end

function state_vector(a::LinkedVorticity3D)
    return state_vector(a.child_object)
end

function update_using_state_vector!(
    a::LinkedVorticity3D,
    state_vect::Vector{Float64})
    return update_using_state_vector!(a.child_object, state_vect)
end

function state_time_derivative(
    a::LinkedVorticity3D,
    inducing_bodies::Vorticity3D)

    return state_time_derivative(a, inducing_bodies)
end

function vorticity_vector_length(a::LinkedVorticity3D)
    return vorticity_vector_length(a.source_object)
end

function vorticity_vector(a::LinkedVorticity3D)
    return vorticity_vector(a.source_object)
end

function update_using_vorticity_vector!(
    a::LinkedVorticity3D,
    vort_vect::Vector{Float64})

    @assert(size(a.vorticity_vector_transform_matrix, 2) ==
        vorticity_vector_length(a.source_object), string("LinkedVorticity3D:",
        " could not update child vorticity vector due to incorrectly shaped ",
        "transformation matrix. Shape of matrix was ", 
        size(a.vorticity_vector_transform_matrix), " where the second index ",
        "should match the length of the sourceobject vorticity vector which ",
        "was ", vorticity_vector_length(a.source_object)))
    child_vort_vec = a.vorticity_vector_transform_matrix * vort_vect
    update_using_vorticity_vector!(a.child_object, child_vort_vec)
    return
end

function vorticity_vector_velocity_influence(
    a::LinkedVorticity3D,
    mes_pnt::Vector3D
    )

    # We don't do the parent's influence!
    inf_mat = vorticity_vector_velocity_influence(a.child_object, mes_pnt)
    trns_mat = a.vorticity_vector_transform_matrix
    @assert(size(inf_mat, 2) == size(trns_mat, 1),
        string("LinkedVorticity3D: The vorticity vector velocity influce",
        " matrix for a single point was of dimensions ", size(inf_mat),
        ". This second index should match the first value of the ",
        "transformation matrix that links the vorticity vector to that of ",
        "the source object. The dimensions of this matrix are ",
        size(trns_mat), "."))
    if any(isnan.(trns_mat)) || any(isinf.(trns_mat))
        @warn "LinkedVorticity3D: Vorticity transformation matrix contained"*
            " NaN or Inf values."
    end
    inf_mat = inf_mat * trns_mat
    return inf_mat
end

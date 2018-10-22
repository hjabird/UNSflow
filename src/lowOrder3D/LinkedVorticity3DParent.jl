#===============================================================================
    LinkedVorticity3DParent.jl

    # DESCRIPTION: TO DO.

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

mutable struct LinkedVorticity3DParent <: Vorticity3DProxy
    source_object :: Vorticity3D
    target_object :: LinkedVorticity3DChild
    vorticity_vector_transform_matrix :: Matrix{Float64}

    function LinkedVorticity3DParent(
        source_object::Vorticity3D, 
        target_object::Vorticity3D, 
        transformation_matrix::Matrix{T}) where T <: Real

        tmat = transformation_matrix
        @assert(size(tmat, 1) == 
            vorticity_vector_length(target_object.source_object),
            string("LinkedVorticity3DParent: the vorticity vector length of the ",
            "child object does not match the number or rows in the ",
            "transformation matrix. Child object vorticity vector length is ",
            vorticity_vector_length(target_object.source_object), 
            " and size of the  matrix is ", size(tmat)))
        @assert(size(tmat, 2) == vorticity_vector_length(source_object),
            string("LinkedVorticity3DParent: the vorticity vector length of the ",
            "source object does not match the number or columns in the ",
            "transformation matrix. Source object vorticity vector length is ",
            vorticity_vector_length(source_object), " and size of the  ",
            "matrix is ", size(tmat)))
        if any(isnan.(tmat)) || any(isinf.(tmat))
            @warn "Transformation matrix contained NaN or Inf values."
        end
        new(source_object, target_object, transformation_matrix)
    end
end


function update_using_vorticity_vector!(
    a::LinkedVorticity3DParent,
    vort_vect::Vector{T}) where T <: Real

    @assert(size(a.vorticity_vector_transform_matrix, 2) ==
        vorticity_vector_length(a.source_object), string("LinkedVorticity3DParent:",
        " could not update child vorticity vector due to incorrectly shaped ",
        "transformation matrix. Shape of matrix was ", 
        size(a.vorticity_vector_transform_matrix), " where the second index ",
        "should match the length of the sourceobject vorticity vector which ",
        "was ", vorticity_vector_length(a.source_object)))
    update_using_vorticity_vector!(a.source_object, vort_vect)
    child_vort_vec = a.vorticity_vector_transform_matrix * vort_vect
    update_using_vorticity_vector!(a.target_object.source_object, child_vort_vec)
    return
end

function vorticity_vector_velocity_influence(
    a::LinkedVorticity3DParent,
    mes_pnt::Vector3D
    )

    inf_mat = vorticity_vector_velocity_influence(
                                        a.target_object.source_object, mes_pnt)
    trns_mat = a.vorticity_vector_transform_matrix
    @assert(size(inf_mat, 2) == size(trns_mat, 1),
        string("LinkedVorticity3DParent: The vorticity vector velocity influce",
        " matrix for a single point was of dimensions ", size(inf_mat),
        ". This second index should match the first value of the ",
        "transformation matrix that links the vorticity vector to that of ",
        "the source object. The dimensions of this matrix are ",
        size(trns_mat), "."))
    if any(isnan.(trns_mat)) || any(isinf.(trns_mat))
        @warn "LinkedVorticity3DParent: Vorticity transformation matrix contained"*
            " NaN or Inf values."
    end
    inf_mat = inf_mat * trns_mat
    inf_mat += vorticity_vector_velocity_influence(a.source_object, mes_pnt)
    return inf_mat
end

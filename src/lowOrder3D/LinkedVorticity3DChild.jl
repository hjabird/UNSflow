#===============================================================================
    LinkedVorticity3DChild.jl

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

mutable struct LinkedVorticity3DChild <: Vorticity3DProxy
    source_object :: Vorticity3D

    function LinkedVorticity3DChild(
        source_object::Vorticity3D)

        new(source_object)
    end
end

function vorticity_vector_length(
    a::LinkedVorticity3DChild)
    # The parent "owns" the controlling vorticity vector.
    return 0
end

function vorticity_vector(a::LinkedVorticity3DChild)
    # The parent "owns" the controlling vorticity vector.
    return zeros(0)
end

function update_using_vorticity_vector!(
    a::LinkedVorticity3DChild,
    vort_vect::Vector{T}) where T <: Real
    @assert(length(vort_vect) == 0, string("LinkedVorticity3DChild:",
        " A LinkedVorticity3DChild's vorticity is set by its parent object. ",
        "The length of its vorticity vector is therefore implicitly zero.",
        " If you really want to do this go update_using_vorticity_vector!",
        "(this.source_object, vort_vector)."))
    return
end

function vorticity_vector_velocity_influence(
    a::LinkedVorticity3DChild,
    mes_pnt::Vector3D
    )

    # Since we want to link the parent and child vorticity indexes, the parent
    # calls the child's source object vorticity_vector_velocity_influence
    # fn. 
    inf_mat = zeros(3, 0)
    return inf_mat
end

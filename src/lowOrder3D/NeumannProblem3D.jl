#===============================================================================
    NeumannProblem3D.jl

    A flow problem with Neumann boundary conditions: ie perscribed velocity
    boundary conditions.

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

mutable struct NeumannProblem3D
    bc_points :: Vector{Vector3D}
    bc_velocities :: Vector{Vector3D}
    variable_vorticities :: Vorticity3D
    invariant_vorticities :: Vorticity3D
end

function solve!(a::NeumannProblem)
    check_valid(a)
    mtrx = normal_velocity_influence_matrix(a.variable_vorticities, 
        a.bc_points, a.bc_velocities)

    fs_vector = map(
        x->dot(induced_velocity(a.invariant_vorticities, x[1]), x[2]), 
        zip(s_centres, s_normals))
        
    # mtrx*solution + fs_vector = 0
    soln = (mtrx) \ (-fs_vector)
    
    update_using_vorticity_vector!(a.variable_aero, soln)
    return
end

function check_valid(a::NeumannProblem)
    len = length(a.bc_points)
    @assert(length(a.bc_velocities) == len,
        string("Length of bc_points and bc_normals vector did not match. ",
            "length(bc_normals) = ", length(a.bc_normals), " and ",
            "length(bc_points) = ", len))
    @assert(all(isfinite.(a.bc_points), "Non-finite points in bc_points."))
    @assert(all(isfinite.(a.bc_velocities)), "Non-finite vectors in "*
            "bc_velocities.")
    @assert(len == vorticity_vector_length(a.variable_vorticities),
        string("Length of boundary condition vectors did not match that ",
        "of the the vorticity control vector: bc vector length was ",
        len, " and vorticity control vector length was ", 
        vorticity_vector_length(a.variable_vorticities)))
    return
end

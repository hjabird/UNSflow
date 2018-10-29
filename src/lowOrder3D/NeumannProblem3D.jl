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
    bc_directions :: Vector{Vector3D}
    bc_velocities :: Vector{Float64}
    variable_vorticities :: Vorticity3D
    invariant_vorticities :: Vorticity3D

    function NeumannProblem3D(
        bc_points :: Vector{Vector3D},
        bc_directions :: Vector{Vector3D},
        bc_velocities :: Vector{Float64},
        variable_vorticities :: Vorticity3D,
        invariant_vorticities :: Vorticity3D)

        new(bc_points, bc_directions, bc_velocities, 
            variable_vorticities, invariant_vorticities)
    end
end

function NeumannProblem3D()
    bc_points = Vector{Vector3D}()
    bc_directions = Vector{Vector3D}()
    bc_velocities = Vector{Float64}()
    variable_vorticities = Vorticity3DSimpleCollector()
    invariant_vorticities = Vorticity3DSimpleCollector()

    NeumannProblem3D(bc_points, bc_directions, bc_velocities, 
        variable_vorticities, invariant_vorticities)
end

function vorticities(a::NeumannProblem3D)
    return Vorticity3DSimpleCollector(
        a.variable_vorticities, a.invariant_vorticities)
end

function solve!(a::NeumannProblem3D)
    check_valid(a)
    mtrx = normal_velocity_influence_matrix(a.variable_vorticities, 
        a.bc_points, a.bc_directions)

    fs_vector = map(
        x->dot(induced_velocity(a.invariant_vorticities, x[1]), x[2]), 
        zip(a.bc_points, unit.(a.bc_directions)))
        
    # mtrx*solution + fs_vector = velocities
    soln = (mtrx) \ (a.bc_velocities - fs_vector)
    
    update_using_vorticity_vector!(a.variable_vorticities, soln)
    return
end

function check_valid(a::NeumannProblem3D)
    len = length(a.bc_points)
    @assert(length(a.bc_directions) == len,
        string("Length of bc_points and bc_directions vector did not match. ",
            "length(bc_directions) = ", length(a.bc_directions), " and ",
            "length(bc_directions) = ", len))
    @assert(length(a.bc_velocities) == len,
        string("Length of bc_points and bc_velocities vector did not match. ",
            "length(bc_velocities) = ", length(a.bc_velocities), " and ",
            "length(bc_velocities) = ", len))
    @assert(all(isfinite.(a.bc_points)), "Non-finite points in bc_points.")
    @assert(all(isfinite.(a.bc_directions)), "Non-finite vectors in "*
            "bc_directions.")
    @assert(all(isfinite.(a.bc_velocities)), "Non-finite values in "*
            "bc_velocities.")
    @assert(len == vorticity_vector_length(a.variable_vorticities),
        string("Length of boundary condition vectors did not match that ",
        "of the the vorticity control vector: bc vector length was ",
        len, " and vorticity control vector length was ", 
        vorticity_vector_length(a.variable_vorticities)))
    return
end

function add_bc!(a::NeumannProblem3D,
    point::Vector3D, direction::Vector3D, velocity::T=0.0) where T <: Real
    if isnan(point) || isinf(point)
        @warn "Point contains NaN or Inf component(s)"
    end
    if isnan(direction) || isinf(direction)
        @warn "direction contains NaN or Inf component(s)"
    end
    if isnan(velocity) || isinf(velocity)
        @warn "velocity is NaN or Inf"
    end
    push!(a.bc_points, point)
    push!(a.bc_directions, diretion)
    push!(a.bc_velocities, velocity)
    return
end

function add_bc!(a::NeumannProblem3D,
    points::Array{Vector3D}, directions::Array{Vector3D}, 
    velocities::Array{T}=fill(0.0, length(points))) where T <: Real

    len = length(points)
    @assert(len == length(directions),
        string("Length of points vector did not match that of directions",
            " vector. length(points) = ", len, " and length(directions) = ",
            length(directions)))
    @assert(len == length(velocities),
        string("Length of points vector did not match that of velocities",
            " vector. length(points) = ", len, " and length(velocities) = ",
            length(velocities)))
    append!(a.bc_points, vec(points))
    append!(a.bc_directions, vec(directions))
    append!(a.bc_velocities, vec(velocities))
    return
end

#===============================================================================
    VortexLatticeMethod.jl

    An implementation of the steady vortex lattice method.

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

mutable struct VortexLatticeMethod
    free_stream :: Vector3D

    variable_aero :: Vorticity3DSimpleCollector
    bc_geometry :: Vector{BilinearQuadSurf}

    function VortexLatticeMethod(
        wing_geometry :: BilinearQuadSurf, 
        free_stream :: Vector3D,
        wake_positions :: Vector{T}) where T <: Real

        wing_aero_original = VortexRingLattice(wing_geometry)
        # The wake aero is er... fun to generate.
        ni, nj = size(wing_geometry)
        wake_points = 
            Matrix{Vector3D}(undef, (length(wake_positions) + 1, nj + 1))
        # Get the points on the TE and put them in every row in our matrix:
        for i = 1 : size(wake_points, 1)
            wake_points[i, :] = wing_geometry.coordinates[end, :]
        end
        # And generate downstream vectors
        wp = vcat(T(0), wake_positions)
        dstream_vects = wp .* [unit(free_stream) for i = 1 : length(wp)]
        for j = 1 : size(wake_points, 2)
            wake_points[:, j] = wake_points[:, j] .+ dstream_vects
        end
        wake_aero = VortexRingLattice(BilinearQuadSurf(wake_points))
        transform_mat = trailing_edge_to_wake_vorticity_transform_matrix(
            wing_aero_original, wake_aero)
        wake_child = LinkedVorticity3DChild(wake_aero)
        aero = Vorticity3DSimpleCollector(
            LinkedVorticity3DParent(wing_aero_original, wake_child, 
                transform_mat),
            wake_child
        )

        new(free_stream, aero, [wing_geometry])
    end
end

function solve!(a::VortexLatticeMethod)
    s_normals = mapreduce(x->vec(normals(x)), vcat, a.bc_geometry)
    s_centres = mapreduce(x->vec(centres(x)), vcat, a.bc_geometry)
    
    @assert(length(s_normals) == vorticity_vector_length(a.variable_aero),
        string("VortexLatticeMethod: ",
            "The number of geometric constraints does not match the ",
            "number of variables in the system. Did you forget to add to ",
            "the boundary conditions of the problem? There are ",
            length(s_normals), " constraint points/normals and ",
            vorticity_vector_length(a.variable_aero), " variables."))

    mtrx = normal_velocity_influence_matrix(a.variable_aero, 
        s_centres, s_normals)

    fs_vector = map(
        x->dot(a.free_stream, x[2]), 
        zip(s_centres, s_normals))
        
    # mtrx*solution + fs_vector = 0
    soln = (mtrx) \ (-fs_vector)
    
    update_using_vorticity_vector!(a.variable_aero, soln)
    return
end

# Make the wake the correct shape.
function relax_wake!(
    a::VortexLatticeMethod; 
    relaxation_factor::T=0.1) where T <: Real

    points = a.variable_aero[2].geometry.coordinates
    for i = 2 : size(points, 1)
        vels = map(
            x-> induced_velocity(a.variable_aero[2], x) + 
                induced_velocity(a.variable_aero[1], x) + a.free_stream,
            points[i, :])
        # Forward difference:
        comparison_vect = unit.(points[i, :] - points[i-1, :])
        # Don't move in the direction of the comparison vector - only normal.
        moves = map(
            x->(x[1] - dot(x[1], x[2]) * x[2]) * relaxation_factor, 
            zip(vels, comparison_vect))
        for j = i : size(points, 1)
            points[j,:] += moves
        end
    end
    return
end

function trailing_edge_to_wake_vorticity_transform_matrix( 
    wing_lattice::VortexRingLattice, wake_lattice::VortexRingLattice)
    # Each strip in the wake has the same ring vorticity. This matches the
    # TE ring the strip is associated with. We want to control the wake the
    # wing vorticity vector - specifically the TE part. We therefore need a 
    # transformation matrix.
    was = wake_lattice.vorticity
    wis = wing_lattice.vorticity

    wing_te_vort_vec_idxs = extract_vorticity_vector_indexes(
        wing_lattice, [size(wis, 1)], collect(1 :size(wis, 2)))
    wake_transform = zeros( length(was), length(wis) )
    for j = 1 : size(was, 2)
        wake_idxs = extract_vorticity_vector_indexes(
            wake_lattice, collect(1:size(was, 1)), [j])
        wake_transform[wake_idxs, wing_te_vort_vec_idxs[j]] .= 1
    end
    return wake_transform
end
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
    wing_geometry :: BilinearQuadSurf
    free_stream :: Vector3D

    wing_aerodynamic :: VortexRingLattice
    wake_aerodynamic :: VortexRingLattice

    function VortexLatticeMethod(
        wing_geometry :: BilinearQuadSurf, 
        free_stream :: Vector3D,
        wake_positions :: Vector{T}) where T <: Real

        wing_aero = VortexRingLattice(wing_geometry)
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
        new(wing_geometry, free_stream, wing_aero, wake_aero)
    end
end

function solve!(a::VortexLatticeMethod)
    wing_normals = vec(normals(a.wing_geometry))
    wing_centres = vec(centres(a.wing_geometry))
    
    wing_mtrx = normal_velocity_influence_matrix(a.wing_aerodynamic, 
        wing_centres, wing_normals)

    wake_mtrx = normal_velocity_influence_matrix(a.wake_aerodynamic,
        wing_centres, wing_normals)
    wake_trnsform = trailing_edge_to_wake_vorticity_transform_matrix(a)
    wake_mtrx *= wake_trnsform

    fs_vector = map(
        x->dot(a.free_stream, x[2]), 
        zip(wing_centres, wing_normals))

    # (wake_mtrx + wing_mtrx)*solution + fs_vector = 0
    soln = (wake_mtrx + wing_mtrx) \ (-fs_vector)
    
    wake_soln = wake_trnsform * soln
    update_using_vorticity_vector!(a.wing_aerodynamic, soln)
    update_using_vorticity_vector!(a.wake_aerodynamic, wake_soln)
    return
end

# Make the wake the correct shape.
function relax_wake!(
    a::VortexLatticeMethod; 
    relaxation_factor::T=0.1) where T <: Real

    points = a.wake_aerodynamic.geometry.coordinates
    for i = 2 : size(a.wake_aerodynamic.geometry.coordinates, 1)
        vels = map(
            x-> induced_velocity(a.wake_aerodynamic, x) + 
                induced_velocity(a.wing_aerodynamic, x) + a.free_stream,
            points[i, :])
        # Forward difference:
        comparison_vect = unit.(points[i, :] - points[i-1, :])
        # Don't move in the direction of the comparison vector - only normal.
        moves = map(
            x->(x[1] - dot(x[1], x[2]) * x[2]) * relaxation_factor, 
            zip(vels, comparison_vect))
        for j = i : size(a.wake_aerodynamic.geometry.coordinates, 1)
            a.wake_aerodynamic.geometry.coordinates[j,:] += moves
        end
    end
    return
end

function trailing_edge_to_wake_vorticity_transform_matrix( 
    a::VortexLatticeMethod)
    # Each strip in the wake has the same ring vorticity. This matches the
    # TE ring the strip is associated with. We want to control the wake the
    # wing vorticity vector - specifically the TE part. We therefore need a 
    # transformation matrix.
    # The vorticity vector will go in the i direction of the wake first.
    wake_veclen = vorticity_vector_length(a.wake_aerodynamic)
    ni, nj = size(a.wing_aerodynamic.geometry)
    nk, nl = size(a.wake_aerodynamic.geometry)
    wake_trnsform = zeros( nk * nl, ni * nj )
    for i = 1 : nj  # Once for each TE ring on Wing
        irange_start = (i - 1) * nk + 1
        irange_end = i * nk
        j_pos = i * ni
        wake_trnsform[irange_start:irange_end, j_pos] = 
            ones(irange_end - irange_start + 1)
    end
    return wake_trnsform
end
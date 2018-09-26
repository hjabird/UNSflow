#===============================================================================
    VortexParticleVolumeAdaptive.jl

    Redistributes particles onto a regular grid within the volume they
    represent.

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
include("Vorticity3DAdaptive.jl")
include("VortexParticle3D.jl")
include("Vortex3DRegularisationFunctions.jl")
include("RedistributionScheme.jl")

mutable struct VortexParticleVolumeAdaptive <: Vorticity3DAdaptive
    children :: Vector{VortexParticle3D}

    gridsize :: Float64
    new_particle_regularisation :: Vortex3DRegularisationFunctions
    redistribution_scheme :: Function

    function VortexParticleVolumeAdaptive(
        particles::Vector{VortexParticle3D},
        gridsize :: Float64,
        new_particle_regularisation=threed_winckelmans_kernels(),
        redistribution_scheme=second_order_redistribution_scheme)

        new(particles, gridsize,
            new_particle_regularisation, redistribution_scheme)
    end
end

function adaptive_update!(a::VortexParticleVolumeAdaptive)
    # Calculate the extrema and how many the number of particles in each dir.
    minmaxx = extrema(map(x->x.geometry.coord.x, a.children))
    minmaxy = extrema(map(x->x.geometry.coord.y, a.children))
    minmaxz = extrema(map(x->x.geometry.coord.z, a.children))
    nx = ceil((minmaxx[2] - minmaxx[1]) / a.gridsize) + 1
    ny = ceil((minmaxy[2] - minmaxy[1]) / a.gridsize) + 1
    nz = ceil((minmaxz[2] - minmaxz[1]) / a.gridsize) + 1
    stepx = (minmaxx[2] - minmaxx[1]) / (nx - 1)
    stepy = (minmaxy[2] - minmaxy[1]) / (ny - 1)
    stepz = (minmaxz[2] - minmaxz[1]) / (nz - 1)
    # Work out the area we need to consider for our redistribution scheme.
    u_crit = mapreduce(x->a.redistribution_scheme(x) > 0 ? 99999 : x,
        min, [0.1:0.1:20;])
    # Box our particles up into u_crit sized boxes
    boxes = Dict{Tuple{Int32, Int32, Int32}, Vector{VortexParticle3D}}()
    function coord_to_box_idx(coord::Vector3D)
        return (ceil(coord.x - minmaxx[1])/u_crit,
                ceil(coord.y - minmaxy[1])/u_crit,
                ceil(coord.z - minmaxz[1])/u_crit)
    end
    for child in a.children
        idx = coord_to_box_idx(child.geometry.coord)
        push!(get!(boxes, idx, Vector{VortexParticle3D}()), child)
    end
    # Now we can build up a new collection of particles.
    new_particles = Vector{VortexParticle3D}()
    function redistributer(coordinate::Vector3D, particle::VortexParticle3D)
        pc = particle.geometry.coord
        mc = coordinate
        rs = a.redistribution_scheme
        coef = rs(abs(pc.x - mc.x)) * rs(abs(pc.y - mc.y))*rs(abs(pc.z - mc.z))
        return coef * particle.vorticity
    end
    function make_new_particle(i::Int64, j::Int64, k::Int64)
        coordinate = Vector3D(minmaxx[1] + stepx * (i - 1),
            minmaxy[1] + stepy * (j - 1), minmaxz[1] + stepz * (k - 1))
        box_idx = coord_to_box_idx(coordinate)
        npidxs = [(i, j, k) for   i = box_idx[1]-1:box_idx[1]+1,
                                j = box_idx[2]-1:box_idx[2]+1,
                                k = box_idx[3]-1:box_idx[3]+1]
        vorticity = Vector3D(0,0,0)
        for npidx in npidxs
            in_the_box =  get(boxes, npidx, Vector{VortexParticle3D}())
            if length(in_the_box) > 0
                vorticity += mapreduce(x->redistributer(coordinate, x),
                    +, in_the_box)
            end
        end
        if abs(vorticity) > 0
            push!(new_particles, VortexParticle3D(coordinate, vorticity,
                1.3 * a.gridsize, a.new_particle_regularisation))
        end
        return
    end
    all_idxs = [(i, j, k) for i=1:nx, j=1:ny, k=1:nz]
    for idx in all_idxs
        make_new_particle(Int64(idx[1]), Int64(idx[2]), Int64(idx[3]))
    end
    a.children = new_particles
    return
end

#- Reimplentation of collection interaction functions ------------------------=#

# We can only push! VortexParticle3Ds.
function Base.push!(a::VortexParticleVolumeAdaptive, b::VortexParticle3D)
    push!(a.children, b)
end

function Base.append!(a::VortexParticleVolumeAdaptive, b)
    @assert(eltype(b) <: VortexParticle3D, "Only VortexParticle3D can be
        added to the VortexParticleVolumeAdaptive.")
    append!(a.children, b)
end

function Base.pop!(a::VortexParticleVolumeAdaptive)
    return pop!(a.children)
end

function Base.isempty(a::VortexParticleVolumeAdaptive)
    return isempty(a.children)
end

function Base.empty!(a::VortexParticleVolumeAdaptive)
    empty!(a.children)
end

#= END VortexFilamentAdaptive ------------------------------------------------=#

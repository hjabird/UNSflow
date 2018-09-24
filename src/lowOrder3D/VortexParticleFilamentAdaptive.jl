#===============================================================================
    VortexParticleFilamentAdaptive.jl

    A spline based vortex filament that supports particle redistribution.

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
import Dierckx

mutable struct VortexParticleFilamentAdaptive <: Vorticity3DAdaptive
    children :: Vector{VortexParticle3D}

    interpolation_order :: Int64
    gridsize :: Float64
    new_particle_regularisation :: Vortex3DRegularisationFunctions
    redistribution_scheme :: Function
end

function adaptive_update(a::VortexParticleFilamentAdaptive)
    # For a spline for each dimension out of our vortex particles.
    x_coords = map(x->x.coord.x, children)
    y_coords = map(x->x.coord.y, children)
    z_coords = map(x->x.coord.z, children)
    index_me = Vector{Float64}(1:length(x_coords))
    spx = Dierckx.Spline1D(index_me, x_coords)
    spy = Dierckx.Spline1D(index_me, y_coords)
    spz = Dierckx.Spline1D(index_me, z_coords)

    # So we need to measure in space along the spline...
    
end

#- Reimplentation of collection interaction functions ------------------------=#

# We can only push! VortexParticle3Ds.
function push!(a::VortexParticleFilamentAdaptive, b::VortexParticle3D)
    push!(a.children, b)
end

function append!(a::VortexParticleFilamentAdaptive, b)
    @assert(eltype(b) <: VortexParticle3D, "Only VortexParticle3D can be
        added to the VortexParticleFilamentAdaptive.")
    append!(a.children, b)
end

function Base.pop!(a::VortexParticleFilamentAdaptive)
    return pop!(a.children)
end

function Base.isempty(a::VortexParticleFilamentAdaptive)
    return isempty(a.children)
end

function Base.empty!(a::VortexParticleFilamentAdaptive)
    empty!(a.children)
end

#= END VortexFilamentAdaptive ------------------------------------------------=#

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
include("RedistributionScheme.jl")
import Dierckx

mutable struct VortexParticleFilamentAdaptive <: Vorticity3DAdaptive
    children :: Vector{VortexParticle3D}

    gridsize :: Float64
    new_particle_regularisation :: Vortex3DRegularisationFunctions
    redistribution_scheme :: Function

    function VortexParticleFilamentAdaptive(
        geometry_definition :: Function, # define geometry in [-1, 1]
        strength :: Float64,
        gridsize :: Float64,
        new_particle_regularisation=threed_winckelmans_kernels(),
        redistribution_scheme=m4_redistribution_scheme)

        @assert(hasmethod(geometry_definition, (Float64)),
            "geometry_definition must be a function that returns the ",
            "coordinate (as a Vector3D) of a point on the spline for ",
            "a Float64 argument in [-1, 1].")
        @assert(0 < gridsize, "gridsize must be positive.")

        dx = x -> ForwardDiff.derivative(geometry_definition, x)


        new(particles, gridsize,
            new_particle_regularisation, redistribution_scheme)
    end
end

function adaptive_update!(a::VortexParticleFilamentAdaptive)
    # For a spline for each dimension out of our vortex particles.
    spline_x = Vector{Float64}(1:length(a.children))
    spx = Dierckx.Spline1D(spline_x, map(x->x.coord.x, a.children))
    spy = Dierckx.Spline1D(spline_x, map(x->x.coord.y, a.children))
    spz = Dierckx.Spline1D(spline_x, map(x->x.coord.z, a.children))
    # Functions for coord and derivative
    xp = x->Vector3D(spx(x), spy(x), spz(x))
    dl = x->Vector3D(derivative(spx, x), derivative(spy, x), derivative(spz, x))
    # Spline for calculating lengths
    gradient_points = Vector{Float64}(1:0.25:length(a.children))
    spl = Dierckx.Spline1D(gradient_points, abs.(dl.(gradient_points)))
    spline_length = integrate(spl, 1, length(a.children))
    # Now we can work on placing the new particles
    new_child_count = ceil(spline_length / a.gridsize)
    new_x_locs = zeros(new_child_count)
    new_x_locs[1] = 1
    for i = 2:new_child_count - 1
        # Newton-Raphson iteration...
        err = 99999
        x = i * (length(x_coords) / length(new_x_locs))
        while err > 0.0001                                     # EVIL CONSTANT!
            x -= (integrate(spl, 0, x) - (i - 1) *
                (spline_length/new_child_count)) / spl(x)
        end
        new_x_locs[i] = x
    end
    new_x_locs[end] = length(a.children)
    # Assign coordinates and size:
    new_children = Vector{VortexParticle3D}()
    for i = 1 : new_child_count
        loc = new_x_locs[i]
        new_children.push!(VortexParticle3D(xp(loc),
            Vector3D(0,0,0), a.gridsize * 0.75), a.new_particle_regularisation)
    end
    # Now we can assign strengths.
    for i = 1:new_child_count
        for j = 1:length(a.children)
            u = abs(integrate(spl, new_x_locs[i], j)) /
                (spline_length / new_child_count)
            fraction = a.redistribution_scheme(u)
            new_children[i].vorticity += a.children[j].vorticity * fraction
        end
    end
    a.children = new_children
    return
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

function funcs_to_equal_spacing(
    coord_fn :: Function,
    coord_deriv_fn :: Function,
    min_arg :: Real,
    max_arg :: Real,
    spacing :: Real
    )

    spline_length = integrate(spl, 1, length(a.children))
    # Now we can work on placing the new particles
    new_child_count = ceil(spline_length / spacing)
    new_x_locs = zeros(new_child_count)
    new_x_locs[1] = min_arg
    for i = 2:new_child_count - 1
        # Newton-Raphson iteration...
        err = 99999
        x = i * (length(x_coords) / length(new_x_locs))
        while err > 0.001 * spacing                             # EVIL CONSTANT!
            x -= (integrate(spl, 0, x) - (i - 1) *
                (spline_length/new_child_count)) / spl(x)
        end
        new_x_locs[i] = x
    end
    new_x_locs[end] = max_arg
end

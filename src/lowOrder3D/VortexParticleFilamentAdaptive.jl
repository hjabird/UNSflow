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
import QuadGK

mutable struct VortexParticleFilamentAdaptive <: Vorticity3DAdaptive
    children :: Vector{VortexParticle3D}

    gridsize :: Float64
    new_particle_regularisation :: Vortex3DRegularisationFunctions
    redistribution_scheme :: Function

    function VortexParticleFilamentAdaptive(
        geometry_definition::Function,
        strength :: Float64,
        gridsize :: Float64,
        new_particle_regularisation=threed_winckelmans_kernels(),
        redistribution_scheme=m4_redistribution_scheme)

        @assert(hasmethod(geometry_definition, Tuple{Float64}),
            "geometry_definition must be a function that returns the ",
            "coordinate (as a Vector3D) of a point on the spline for ",
            "a Float64 argument in [-1, 1].")
        @assert(0 < gridsize, "gridsize must be positive.")

        spline_x, spline_y, spline_z = three_d_function_to_splines(
            geometry_definition, -1., 1., 0.02)
        dx = x -> Vector3D(
            Dierckx.derivative(spline_x, x),
            Dierckx.derivative(spline_y, x),
            Dierckx.derivative(spline_z, x))
        new_x_locs, spline_length =
            funcs_to_equal_spacing(dx, -1., 1., gridsize)
        particles = Vector{VortexParticle3D}()
        actual_spacing = spline_length / (length(new_x_locs) - 1)
        for i = 1 : length(new_x_locs)
            p = VortexParticle3D(geometry_definition(new_x_locs[i]),
                unit(dx(new_x_locs[i])) * strength * actual_spacing,
                actual_spacing * 1.3)
            push!(particles, p)
        end
        new(particles, gridsize,
            new_particle_regularisation, redistribution_scheme)
    end
end

function adaptive_update!(a::VortexParticleFilamentAdaptive)
    # For a spline for each dimension out of our vortex particles.
    println("In fn.")
    spline_x = Vector{Float64}(1:length(a.children))
    spx = Dierckx.Spline1D(spline_x, map(x->x.geometry.coord.x, a.children))
    spy = Dierckx.Spline1D(spline_x, map(x->x.geometry.coord.y, a.children))
    spz = Dierckx.Spline1D(spline_x, map(x->x.geometry.coord.z, a.children))
    # Functions for coord and derivative
    xp = x->Vector3D(spx(x), spy(x), spz(x))
    dl = x->Vector3D(Dierckx.derivative(spx, x), Dierckx.derivative(spy, x),
        Dierckx.derivative(spz, x))
    # Spline for calculating lengths
    println("Point1")
    gradient_points = Vector{Float64}(1:0.25:length(a.children))
    spl = Dierckx.Spline1D(gradient_points, abs.(dl.(gradient_points)))
    println("Point2")
    new_x_locs, spline_length =
                funcs_to_equal_spacing(dl, 1, length(a.children), a.gridsize)
    println("Point3")
    # Assign coordinates and size:
    new_children = Vector{VortexParticle3D}()
    for i = 1 : length(new_x_locs)
        loc = new_x_locs[i]
        push!(new_children, VortexParticle3D(xp(loc),
            Vector3D(0,0,0), a.gridsize * 0.75, a.new_particle_regularisation))
    end
    # Now we can assign strengths.
    actual_spacing = spline_length / (length(new_x_locs) - 1)
    for i = 1 : length(new_x_locs)
        for j = 1:length(a.children)
            u = abs(Dierckx.integrate(spl, new_x_locs[i], j)) /  actual_spacing
            fraction = a.redistribution_scheme(u)
            new_children[i].vorticity += a.children[j].vorticity * fraction
        end
    end
    a.children = new_children
    return
end

#- Reimplentation of collection interaction functions ------------------------=#

# We can only push! VortexParticle3Ds.
function Base.push!(a::VortexParticleFilamentAdaptive, b::VortexParticle3D)
    push!(a.children, b)
end

function Base.append!(a::VortexParticleFilamentAdaptive, b)
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
    coord_deriv_fn :: Function,
    min_arg :: Real,
    max_arg :: Real,
    spacing :: Real
    )
    dlength = x::Float64 -> abs(coord_deriv_fn(x))
    spline_length, _ = QuadGK.quadgk(dlength, min_arg, max_arg)
    # Now we can work on placing the new particles
    new_child_count = ceil(spline_length / spacing) + 1
    actual_spacing = spline_length / (new_child_count - 1)
    new_x_locs = zeros(Int32(new_child_count))
    new_x_locs[1] = min_arg
    println("Spacing 2")
    for i = 2:new_child_count - 1
        # Newton-Raphson iteration...
        err = 99999
        x = i * ((max_arg-min_arg) / length(new_x_locs)) - min_arg
        iter_count = 0
        while err > 0.001 * spacing                             # EVIL CONSTANT!
            x_minus1 = x
            println("Prequad")
            x -= (QuadGK.quadgk(dlength, min_arg, x)[1] - (i - 1) *
                actual_spacing) / dlength(x)
            println("Postquad")
            err = abs(x - x_minus1)
        end
        new_x_locs[Int64(i)] = x
        iter_count += 1
        if iter_count > 30
            println("Iteration limit reached: bombing out with ", err / spacing,
            "error.")
        end
    end
    new_x_locs[end] = max_arg
    println("Spacing end")
    return new_x_locs, spline_length
end

function three_d_function_to_splines(func::Function, min_arg::Real,
        max_arg::Real, resolution :: Real)
    is = Vector{Float64}(min_arg:resolution:max_arg)
    vals = func.(is)
    xs = map(x->x.x, vals)
    ys = map(x->x.y, vals)
    zs = map(x->x.z, vals)
    spx = Dierckx.Spline1D(is, xs)
    spy = Dierckx.Spline1D(is, ys)
    spz = Dierckx.Spline1D(is, zs)
    return spx, spy, spz
end

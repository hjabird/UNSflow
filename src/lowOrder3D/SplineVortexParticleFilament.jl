#===============================================================================
    SplineVortexParticleFilament.jl

    A spline based vortex filament intended to as a means to generate
    vortex particles.

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
import Dierckx

mutable struct SplineVortexParticleFilament <: Vorticity3DCollector
    children :: Vector{VortexParticle3D}
    vorticity :: Float64

    # SPLINES in [-1, 1]
    function SplineVortexParticleFilament(
        xs :: Dierckx.Spline1D, ys :: Dierckx.Spline1D,
        zs :: Dierckx.Spline1D, particle_radius :: Float64,
        strength :: Float64 = 1.0,
        kernel :: Vortex3DRegularisationFunctions =
            threed_winckelmans_kernels()
        )
        @assert(strength != 0.0)
        @assert(particle_radius >= 0)

        x_p = x -> Vector3D(xs(x), ys(x), zs(x))
        dx_p = x -> Vector3D(
            Dierckx.derivative(xs, x),
            Dierckx.derivative(ys, x),
            Dierckx.derivative(zs, x))
        absdx_p = x -> abs(dx_p(x))

        stepping = 0.05;
        eval_locs = collect(-1.:stepping:1.)
        abs_dvals = absdx_p.(eval_locs)
        spl_len = (sum(abs_dvals[1:end-1]) + sum(abs_dvals[2:end])) .*
            stepping .* 0.5
        n_particles = ceil(spl_len / (2 * particle_radius))
        particle_sep = spl_len / n_particles

        particles = Vector{VortexParticle3D}(undef, n_particles)
        for i = 1 : length(particles)
            # Work out the location.
            length_pos = (i - 0.5) * particle_sep
            j = 0
            while stepping * (sum(abs_dvals[1:j]) - 0.5  * abs_dvals[1]) < length_pos
                j += 1
            end
            px = j - (sum(abs_dvals[1:j]) - length_pos) / abs_dvals[j]
            pos = x_p(px)
            # Vorticity is easier.
            spline_dir = dx_p(px)
            vort = particle_sep * strength / abs(dx_p) * dx_p
            particles[i] = VortexParticle3D(pos, vort, stepping/2, kernel)
        end
        new(particles, strength)
    end
end

function centre(this::SplineVortexParticleFilament)
    return centre(VortexParticleSimpleCollector(this.children))
end

function effective_radius(this::SplineVortexParticleFilament)
    return effective_radius(VortexParticleSimpleCollector(this.children))
end

function vorticity(this::SplineVortexParticleFilament)
    return vorticity(VortexParticleSimpleCollector(this.children))
end

function euler!(
    this::SplineVortexParticleFilament,
    influence_field::Vorticity3D,
    dt::Real)
    error("SplineVortexParticleFilament is intended only for the "*
        "intitialisation"*
        " of vortex particles on a given geometry with a single vorticity "*
        "handel. Consider expanding this collector into a "*
        "Vorticity3DSimpleCollector or, if you wanted something that would "*
        "redistribute vortex particles later or, "*
        "VortexParticleFilamentAdaptive.")
    return
end

function state_vector_length(a::SplineVortexParticleFilament)
    error("SplineVortexParticleFilament is intended only for the "*
        "intitialisation"*
        " of vortex particles on a given geometry with a single vorticity "*
        "handel. Consider expanding this collector into a "*
        "Vorticity3DSimpleCollector or, if you wanted something that would "*
        "redistribute vortex particles later or, "*
        "VortexParticleFilamentAdaptive.")
    return
end

function state_vector(a::SplineVortexParticleFilament)
    error("SplineVortexParticleFilament is intended only for the "*
        "intitialisation"*
        " of vortex particles on a given geometry with a single vorticity "*
        "handel. Consider expanding this collector into a "*
        "Vorticity3DSimpleCollector or, if you wanted something that would "*
        "redistribute vortex particles later or, "*
        "VortexParticleFilamentAdaptive.")
    return
end

function update_using_state_vector!(
    this::SplineVortexParticleFilament,
    state_vect::Vector{Float64})
    error("SplineVortexParticleFilament is intended only for the "*
        "intitialisation"*
        " of vortex particles on a given geometry with a single vorticity "*
        "handel. Consider expanding this collector into a "*
        "Vorticity3DSimpleCollector or, if you wanted something that would "*
        "redistribute vortex particles later or, "*
        "VortexParticleFilamentAdaptive.")
    return
end

function state_time_derivative(
    this::SplineVortexParticleFilament,
    inducing_bodies::Vorticity3D)
    error("SplineVortexParticleFilament is intended only for the "*
        "intitialisation"*
        " of vortex particles on a given geometry with a single vorticity "*
        "handel. Consider expanding this collector into a "*
        "Vorticity3DSimpleCollector or, if you wanted something that would "*
        "redistribute vortex particles later or, "*
        "VortexParticleFilamentAdaptive.")
    return
end

function vorticity_vector_length(this::SplineVortexParticleFilament)
    return 1
end

function vorticity_vector(this::SplineVortexParticleFilament)
    return [this.vorticity]
end

function update_using_vorticity_vector!(
    this::SplineVortexParticleFilament,
    vort_vect::Vector{Float64})
    old_v = this.vorticity
    if old_v != 0
        ratio = vort_vect[1] / old_v
        for child in this.children
            child.vorticity *= ratio
        end
    else
        error("Vorticity was zero. SplineVortexParticleFilament.vorticity "*
            "should never be set to zero to avoid division by zeros later on.")
    end
    return
end

function vorticity_vector_velocity_influence(
    this::SplineVortexParticleFilament,
    mes_pnt::Vector3D
    )
    v = zeros(3, 1)
    for child in this.children
        #= We need to correct for the direction of the vorticity and
        for the relative strength of the filament and the particle it generated
        =#
        abs_str = abs(child.vorticity) / this.vorticity
        v_p = vorticity_vector_velocity_influence(child, mes_pnt) .*
            child.vorticity * abs_str
        v += v_p
    end
    return v
end

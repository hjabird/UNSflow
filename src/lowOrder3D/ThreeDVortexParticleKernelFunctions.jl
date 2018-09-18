#===============================================================================
    ThreeDVortexParticleKernelFunctions.jl

    These methods return the reduction_factor_fn and vorticity_fraction_fn
    used by ThreeDVortexParticle in velocity etc. calculations.

    These are taken from:

    Robertson, Joo and Reich: "Vortex Particle
    Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
    51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
    Materials Conference

    Winckelmans et al,"Vortex methods and their application to trailing wake
    vortex simulations", Comptes Rendus Physique, 2005

    The Winckelmans nomenclature has been adopted for the g(rho), G(rho)
    and zeta(rho).

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

type ThreeDVortexParticleKernelFunctions
    g :: Function
    G :: Function
    zeta :: Function
    radius_modifier :: Function
end

# THE INTERFACE:

# g(rho::Real)
#   The g function defined by Winckelmans. returns a real.

# G(rho::Real)
#   The g function defined by Winckelmans. returns a real.

# zeta(rho::Real)
#   The g function defined by Winckelmans. returns a real.

# radius_modifier()
#   returns a real representing the radius of the vortex particle adjusted
#   for the fact that the support of the regularisation function may not be
#   compact.


function threed_particle_singularity_kernels()
    function vortex_g_singularity(rho::Real)
        return 1.0
    end
    function vortex_G_singularity(rho::Real)
        error("G is not well defined for a singularity.")
    end
    function vortex_zeta_singularity(rho::Real)
        return 0.0
    end
    function radius_modifier()
        return 1.0
    end
    return ThreeDVortexParticleKernelFunctions(
        vortex_g_singularity,
        vortex_G_singularity,
        vortex_zeta_singularity,
        radius_modifier
    )
end


function threed_planetary_kernels()
    function vortex_g_planetary(rho::Real)
        return (min(1, rho)) ^ 3
    end
    function vortex_G_singularity(rho::Real)
        error("Not yet implemented. Please consider contributing this.")
    end
    function vortex_zeta_planetary(rho::Real)
        return rho <  1 ? 3 : 0
    end
    function radius_modifier()
        return 1.0
    end
    return ThreeDVortexParticleKernelFunctions(
        vortex_g_planetary,
        vortex_G_planetary,
        vortex_zeta_planetary,
        radius_modifier
    )
end


function threed_exponential_kernels()
    function vortex_g_exponential(rho::Real)
        return 1 - vortex_f_exponential(rho) / 3.0
    end
    function vortex_G_exponential(rho::Real)
        error("Not yet implemented. Please consider contributing this.")
    end
    function vortex_zeta_exponential(rho::Real)
        return 3.0 * exp(- rho ^ 3)
    end
    function radius_modifier()
        error("Not yet implemented. Please consider contributing this.")
    end
    return ThreeDVortexParticleKernelFunctions(
        vortex_g_exponential,
        vortex_G_exponential,
        vortex_zeta_exponential,
        radius_modifier
    )
end


function threed_winckelmans_kernels()
    function vortex_g_winckelmans(rho::Real)
        a = rho^2 + 2.5
        b = rho^3
        c = rho^2 + 1
        d = a * b
        e = c ^(5./2.)
        return d / e
    end
    function vortex_G_winckelmans(rho::Real)
        return 0.5 * (log(rho^2 + 1) + (rho^2) / (rho^2 + 1))
    end
    function vortex_f_winckelmans(rho::Real)
        a = 15. / 2.
        b = rho ^ 2 + 1
        return a / (b ^ (7./2.))
    end
    function radius_modifier()
        error("Not yet implemented. Please consider contributing this.")
    end
    return ThreeDVortexParticleKernelFunctions(
        vortex_g_winckelmans,
        vortex_G_winckelmans,
        vortex_zeta_winckelmans,
        radius_modifier
    )
end


function threed_tanh_kernels()
    function vortex_g_tanh(rho::Real)
        return tanh(rho^3)
    end
    function vortex_G_tanh(rho::Real)
        error("Not yet implemented. Please consider contributing this.")
    end
    function vortex_zeta_tanh(rho::Real)
        return 3 * sech(rho^3)^2
    end
    function radius_modifier()
        error("Not yet implemented. Please consider contributing this.")
    end
    return ThreeDVortexParticleKernelFunctions(
        vortex_g_tanh,
        vortex_G_tanh,
        vortex_zeta_tanh,
        radius_modifier
    )
end


function threed_gaussian_kernels()
    function vortex_g_gaussian(rho::Float64)
        return erf(rho / sqrt(2.)) - rho * vortex_f_gaussian(rho)
    end
    function vortex_G_gaussian(rho::Real)
        return 0.5 * (log(rho^2/2) + erf(rho^2/2))
    end
    function vortex_zeta_gaussian(rho::Float64)
        return sqrt(2 / pi) * exp((- rho ^2) / 2)
    end
    function radius_modifier()
        error("Not yet implemented. Please consider contributing this.")
    end
    return ThreeDVortexParticleKernelFunctions(
        vortex_g_gaussian,
        vortex_G_gaussian,
        vortex_zeta_gaussian,
        radius_modifier
    )
end


function threed_super_gaussian_kernels()
    function vortex_g_super_gaussian(rho::Float64)
        return erf(rho / sqrt(2.)) - ((2 - rho^2) / (5 - rho^2)) *
            vortex_f_super_gaussian(rho)
    end
    function vortex_G_super_gaussian(rho::Real)
        return 0.5 * (log(rho^2/2) + erf(rho^2/2) - exp(- rho^2 / 2))
    end
    function vortex_zeta_super_gaussian(rho::Float64)
        return sqrt(2 / pi) * (2.5 - rho^2 / 2) * exp((- rho ^2) / 2)
    end
    function radius_modifier()
        error("Not yet implemented. Please consider contributing this.")
    end
    return ThreeDVortexParticleKernelFunctions(
        vortex_g_super_gaussian,
        vortex_G_super_gaussian,
        vortex_zeta_super_gaussian,
        radius_modifier
    )
end
# END ThreeDVortexParticleKernelFunctions =====================================#

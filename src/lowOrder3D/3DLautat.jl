function VortexParticleWakeLAUTATSolution(
    starting_solution :: VortexParticleWakeLAUTATSolution,
    kinematics :: KinemDef3D,
    nsteps :: Int64,
    dt :: Float64
    )

    # First we need to set up methods by which the vortex particle wake can be
    # affected by the wing section. Our discretisation is two straight
    # filaments.

    for i = 1 : nsteps
        euler_forward_step!(starting_solution.wake, dt)
        add_new_particles(starting_solution)
        compute_fourier_terms(starting_solution)
    end
    return starting_solution
end

"""Covert x in [-1, 1] to theta in [0, pi]"""
function x_to_theta(x :: Float64)
    @assert(abs(x) <= 1.0)
    return acos(-x)
end

"""Convert theta in [0, pi] to x in [-1, 1]"""
function theta_to_x(theta :: Float64)
    @assert(0 <= theta <= pi)
    return -cos(theta)
end

""" Obtain geometric information about a point on the wing.

xs and strip_idxs are matched by index. Returns vectors x, chord_direction,
no-camber wing normal direction and detadx times the c direction.
"""
function xs_and_strip_idxs_to_point_surface_info(
    wing :: StripDefinedWing,
    xs :: Vector{Float64},
    strip_idxs :: Vector{Int}
    )

    @assert(length(xs) == length(strip_idxs))
    @assert(all(abs.(xs) .<= 1.0))
    @assert(all(0.0 .< strip_idxs .<= length(wing.strips())))
    n = length(xs)
    detadx_fn = get_surface_detadx_dot_c_fn(wing)
    chord_dir_fn = get_chord_dir_fn(wing)
    surface_fn = get_surface_fn(wing)
    snx, sny, snz = nocamber_normal_splines(wing)
    normal_dir_fn = s -> ThreeDVector(snx(s), sny(s), snz(s))
    x = map(x->surface_fn(x[1], x[2]), zip(strip_idxs, xs))
    c_dir = map(x->chord_dir_fn(x[1], x[2]), zip(strip_idxs, xs))
    n_dir = map(x->normal_dir_fn(x[1]), strip_idxs)
    detadx_c = map(x->detadx_fn(x[1], x[2]), zip(strip_idxs, xs))
    return x, c_dir, n_dir, detadx_c
end


"""Obtain a single index from a strip index and fourier term number."""
function fourier_term_linear_index(
    sim :: VortexParticleWakeLAUTATSolution,
    strip_idx :: Int64,
    fourier_n :: Int64
    )
    @assert(0 <= fourier_n < sim.n_fourier_terms)
    @assert(1 <= strip_idx <= length(sim.wing.strips))
    return (strip_idx - 1) * sim.n_fourier_terms + fourier_n + 1
end

"""Obtain strip index and fourier term number respectively from a linear index.
"""
function fourier_term_vecvec_index(
    sim :: VortexParticleWakeLAUTATSolution,
    linear_idx :: Int64
    )
    @assert(1 <= linear_idx <= sim.n_fourier_terms * length(sim.wing.strips))
    const nf = sim.n_fourier_terms
    f = linear_idx % nf - 1
    s = (linear_idx - f - 1) / nf
    @assert(Float64(s) == Int(s))
    return s, f
end

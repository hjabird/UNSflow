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
    return fourier_term_linear_index(sim.n_fourier_terms,
        length(sim.wing.strips))
end


function fourier_term_linear_index(
    max_n :: Int64,
    n_strips :: Int64,
    fourier_n :: Int64,
    strip_idx :: Int64
    )
    @assert(0 <= fourier_n <= n_fourier_terms)
    @assert(1 <= strip_idx <= n_strips)
    ret = (strip_idx - 1) * (max_n + 1) + fourier_n + 1
    return ret
end


"""Obtain strip index and fourier term number respectively from a linear index.
"""
function fourier_term_vecvec_index(
    sim :: VortexParticleWakeLAUTATSolution,
    linear_idx :: Int64
    )
    return fourier_term_vecvec_index(sim.n_fourier_terms,
        length(sim.wing.strips))
end


function fourier_term_vecvec_index(
    n_fourier_terms :: Int64, n_strips :: Int64,
    linear_idx :: Int64
    )
    @assert(1 <= linear_idx <= n_fourier_terms * n_strips)
    const nf = n_fourier_terms
    f = linear_idx % nf - 1
    s = (linear_idx - f - 1) / nf
    @assert(Float64(s) == Int(s))
    return s, f
end


"""
Returns a function describing the vorticity fourier vorticity density function
with repect to nth term. In terms of theta in [0, pi]. For Coeff = 1.0
"""
function vorticity_density_theta_fn(n :: Int64)
    @assert(n >= 0)
    if n > 0
        f = function(x ::Float64)
            @assert(0 <= x <= Float64(pi))
            r = sin(n * x)
            @assert(isfinite(r))
            return r
        end
    else
        f = function(x :: Float64)
            @assert(0 < x <= Float64(pi))
            r = (1 + cos(x)) / sin(x)
            @assert(isfinite(r))
            return r
        end
    end
    s = x::Float64->f(x)
    @assert(isfinite(s(1.57)))
    return s
end


"""
Returns a function describing the vorticity fourier vorticity density function
with repect to nth term. In terms of x in [-1, 1]. For Coeff = 1.0
"""
function vorticity_density_x_fn(n :: Int64)
    f = vorticity_density_theta_fn(n)
    return x::Float64->f(x_to_theta(x))
end


function vorticity_density_theta_derivative_fn(n :: Int64)
    @assert(n >= 0)
    if n > 0
        f = function(x ::Float64)
            @assert(0 <= x <= Float64(pi))
            r = n * cos(n * x)
            @assert(isfinite(r))
            return r
        end
    else
        f = function(x :: Float64)
            @assert(0 < x <= Float64(pi))
            r = -cos(x) / (sin(x)^2)
            @assert(isfinite(r))
            return r
        end
    end
    s = x::Float64->f(x)
    @assert(isfinite(s(1.57)))
    return s
end


function vorticity_density_x_derivative_fn(n :: Int64)
    f = vorticity_density_theta_fn(n)
    return x::Float64->f(x_to_theta(x)) / sqrt(1 - x^2)
end


function vorticity_density_x(fourier_coeffs :: Vector{Float64}, x :: Float64)
    @assert(-1 < x <= 1)
    @assert(all(isfinite.(fourier_coeffs)))
    a = mapreduce(i->vorticity_density_x_fn(i) * fourier_coeffs[i+1],
        +, 0.0, 0 : length(fourier_coeffs) - 1)
    @assert(isfinite(a))
    return a
end


function vorticity_density_x_derivative(
    fourier_coeffs :: Vector{Float64},
    x :: Float64)
    @assert(-1 < x <= 1)
    @assert(all(isfinite.(fourier_coeffs)))
    a = mapreduce(i->vorticity_density_x_derivative_fn(Int64(i))(x) * fourier_coeffs[i+1],
        +, 0.0, 0 : length(fourier_coeffs) - 1)
    @assert(isfinite(a))
    return a
end


"""
Zeros the wing vorticity and sets a single strip to a fourier function.
"""
function set_wing_to_single_bv_fn!(
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    wing :: StripDefinedWing,
    filament_locs :: Vector{Float64},
    n:: Int64,
    strip_idx :: Int64)

    @assert(n >= 0)
    @assert(0 < strip_idx <= length(fil_wing.filaments_ym))
    @assert(all(-1 .< filament_locs .< 1))
    @assert(length(wing.strips) == length(fil_wing.filaments_yp))

    zero_vorticities!(fil_wing)
    fn = vorticity_density_x_fn(n)
    @assert(isfinite(fn(0.0)))
    add_vorticity!(wing, filament_locs, fil_wing, strip_idx, fn)
    return
end

"""
Set wing vorticity according to coeffs
"""
function set_wing_to_fourier_set!(
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    wing :: StripDefinedWing,
    filament_locs :: Vector{Float64},
    A_vals :: Vector{Vector{Float64}}
    )
    @assert(length(A_vals) == length(fil_wing.filaments_ym))
    @assert(length(fil_wing.filaments_yp)== length(fil_wing.filaments_ym))
    zero_vorticities!(fil_wing)
    for strip = 1 : length(A_vals)
        for n = 0 : length(A_vals[strip]) - 1
            fn = x->vorticity_density_x_fn(n)(x) * A_vals[strip][n+1]
            add_vorticity!(wing, filament_locs, fil_wing, strip, fn)
        end
    end
    return
end

"""
Obtain the wash function at a point due to some induced velocity
"""
function downwash_2D(
    normal_dir_fn :: Function,
    deta_dx_dot_c_fn :: Function,
    surf_fn :: Function,
    ind_vel_fn :: Function,
    spanwise_pos :: Float64,
    chordwise_pos :: Float64
    )

    @assert(abs(chordwise_pos) <= 1)
    loc = surf_fn(spanwise_pos, chordwise_pos)
    cdir = normal_dir_fn(spanwise_pos)
    deta = deta_dx_dot_c_fn(spanwise_pos, chordwise_pos)
    v = ind_vel_fn(loc)
    return dot(v, cdir + deta)
end


"""
Computes a set of fourier integrals over a given strip given a ind_vel_fn.
"""
function downwash_2D_fourier_integrals(
    normal_dir_fn :: Function,
    deta_dx_dot_c_fn :: Function,
    surf_fn :: Function,
    ind_vel_fn :: Function,
    spanwise_pos :: Int64,
    n_max :: Int64,
    ref_vel :: Float64,
    filament_locations :: Vector{Float64}
    )
    @assert(n_max >= 0)
    @assert(ref_vel != 0)

    funcs = get_fourier_integrand_vector(n_max, ref_vel)
    integrals = Vector{Float64}(n_max + 1)
    theta = x_to_theta.(filament_locations)
    wash = map(t->downwash_2D(normal_dir_fn, deta_dx_dot_c_fn, surf_fn,
                ind_vel_fn, Float64(spanwise_pos), theta_to_x(t)), theta)
    for i = 0 : n_max
        integrand = map(x->funcs[i+1](x[1],x[2]), zip(theta, wash))
        s = Spline1D(theta, integrand, k=2)
        integrals[i+1] = integrate(s, 0, pi)
    end
    if !all(isfinite.(integrals))
        error(string("Non-finite integrals! Not good. Computing for span pos ",
            "= ", spanwise_pos, " and n_max = ", n_max,
            ". Integrals: ", integrals))
    end
    return integrals
end


"""
Get integrands for computing the fourier terms. These take theta
and downwash as arguments.
"""
function get_fourier_integrand_vector(n_max :: Int64, ref_vel :: Float64)
    funcs = Vector{Function}(n_max + 1)
    let
        mult = -1. / (ref_vel * pi)
        function kernel(theta :: Float64, wash :: Float64)
            @assert(0. <= theta <= pi)
            return mult * wash
        end
        funcs[1] = kernel
    end
    for i = 1 : n_max
        mult = 2. / (ref_vel * pi)
        function kernel(theta :: Float64, wash :: Float64)
            @assert(0. <= theta <= pi)
            return mult * wash * cos(i * theta)
        end
        funcs[i + 1] = kernel
    end
    return funcs
end


"""
Computes fourier integrals for all the strips given an ind_vel_fn
"""
function downwash_2D_fourier_integrals_all_strips(
    normal_dir_fn :: Function,
    deta_dx_dot_c_fn :: Function,
    surf_fn :: Function,
    ind_vel_fn :: Function,
    spanwise_pos_max :: Int64,
    n_max :: Int64,
    ref_vel :: Float64,
    filament_locs :: Vector{Float64}
    )
    integrals = Vector{Float64}((n_max + 1) * spanwise_pos_max)

    for i = 1 : spanwise_pos_max
        idxs =
            fourier_term_linear_index.(n_max, spanwise_pos_max, 0:n_max, i)
        integrals[idxs] =
            downwash_2D_fourier_integrals(
                normal_dir_fn, deta_dx_dot_c_fn, surf_fn,
                ind_vel_fn, i, n_max, ref_vel, filament_locs
                )
    end
    return integrals
end


"""
Compute self_influence_matrix
"""
function downwash_2D_fourier_integrals_all_strips(
    normal_dir_fn :: Function,
    deta_dx_dot_c_fn :: Function,
    surf_fn :: Function,
    spanwise_pos_max :: Int64,
    n_max :: Int64,
    ref_vel :: Float64,
    setup_fn :: Function,
    filament_locs :: Vector{Float64}
    )

    nt = (n_max + 1) * spanwise_pos_max
    integral_matrix = zeros(nt, nt)

    for s = 1 : spanwise_pos_max
        for n = 0 : n_max
            ind_vel_fn = setup_fn(s, n)
            lindex = fourier_term_linear_index(
                n_max, spanwise_pos_max, n, s)
            integral_matrix[:, lindex] =
                downwash_2D_fourier_integrals_all_strips(
                    normal_dir_fn, deta_dx_dot_c_fn, surf_fn,
                    ind_vel_fn, spanwise_pos_max, n_max,
                    ref_vel, filament_locs
                )
        end
    end

    return integral_matrix
end


function downwash_2D_pure_2D_self_fourier_integrals_matrix(
    chord :: WingChordSection,
    filament_positions :: Vector{Float64},
    n_max :: Int64,
    ref_vel :: Float64
    )
    # Source and measurement locations
    sx = filament_positions
    sy = chord.camber_line(filament_positions)
    mx = filament_positions
    mtheta = x_to_theta.(mx)
    my = chord.camber_line.(mx)
    detadx = x->derivative(chord.camber_line, x)

    u = (dx, dy)-> (dx<=eps(Float64) && dy<=eps(Float64) ? 0 : dy / (2 * pi * (dy^2 + dx^2)))
    v = (dx, dy)-> (dx<=eps(Float64) && dy<=eps(Float64) ? 0 :- dx / (2 * pi * (dy^2 + dx^2)))
    inf = (x, y, xm, ym)-> u(xm-x, ym-y) * detadx(xm) - v(xm-x, ym-y)
    inf_mat = zeros(length(mx), length(sx))
    for i = 1 : length(mx)
        inf_mat[i, :] = map(x->inf(x[1],x[2],mx[i],my[i]), zip(sx, sy))
    end
    if !all(isfinite.(inf_mat))
        error("Non-finite values in pure 2D wash influence matrix!")
    end

    funcs = get_fourier_integrand_vector(n_max, ref_vel)
    a_n = zeros(n_max + 1, n_max + 1)
    for s = 0 : n_max
        density_func_s = vorticity_density_theta_fn(s)
        for m = 0 : n_max
            density_func = vorticity_density_x_fn(m)
            vort = lump_vorticities(chord, density_func, filament_positions)
            kernal = inf_mat * vort .* density_func_s.(mtheta)
            ws = Spline1D(mtheta, kernal, k=2)
            a_n[s + 1, m + 1] = integrate(ws, 0, pi)
        end
    end
    return a_n
end


function downwash_2D_pure_2D_self_fourier_integrals_matrix(
    wing :: StripDefinedWing,
    filament_positions :: Vector{Float64},
    n_max :: Int64,
    ref_vel :: Float64
    )
    dw = zeros((n_max+1)*length(wing.strips), (n_max+1)*length(wing.strips))
    for s = 1 : length(wing.strips)
        idxs = fourier_term_linear_index.(
            n_max, length(wing.strips), 0:n_max, s)
        dw[idxs, idxs] = downwash_2D_pure_2D_self_fourier_integrals_matrix(
            wing.strips[s], filament_positions, n_max,
            ref_vel )
    end
    return dw
end


function wing_velocity_to_fourier_integrals_vector(
    normal_dir_fn :: Function,
    deta_dx_dot_c_fn :: Function,
    untransformed_wing :: StripDefinedWing,
    kinematics :: ThreeDCoordinateTransform,
    spanwise_pos_max :: Int64,
    n_max :: Int64,
    ref_vel :: Float64,
    filament_locations :: Vector{Float64}
    )
    integrals = zeros((n_max + 1) * spanwise_pos_max)
    surf = get_surface_fn(untransformed_wing)
    ind_vel = x->-derivative2(kinematics, x)

    for s = 1 : spanwise_pos_max
        idxs =
            fourier_term_linear_index.(n_max, spanwise_pos_max, 0:n_max, s)
        integrals[idxs] =
            downwash_2D_fourier_integrals(
                normal_dir_fn, deta_dx_dot_c_fn, surf,
                ind_vel, s, n_max, ref_vel, filament_locations
                )
    end

    return integrals
end


function compute_new_particles_and_wing_fourier_integral_matrix(
    wing :: StripDefinedWing,
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    fil_locs :: Vector{Float64},
    old_bound_vorticities :: Vector{Float64},
    ref_vel :: Float64,
    new_particles :: ThreeDVortexParticleSet,
    kelvin_particle_indexes :: Vector{Int64},
    dt :: Float64,
    ind_vel_external :: Function,
    n_max :: Int64
    )
    @assert(length(wing.strips) > 0)
    @assert(length(wing.strips) == length(fil_wing.filaments_yp))
    @assert(isfinite(ref_vel))
    @assert(length(new_particles) ==
                length(kelvin_particle_indexes) + length(wing.strips) * 2 + 1)
    @assert(all(-1 .< fil_locs .< 1))
    @assert(n_max > 0)
    @assert(dt > 0)

    ind_vel_np :: Function = x->x # Induced vel excluding new particles. [placeholder for scope]
    ind_vel_pw :: Function = x->ind_vel(fil_wing, x) + ind_vel(particles, x) # Induced vel due to new particles and wing. [placeholder for scope]
    function reset(sti :: Int64, n :: Int64)
        @assert(1 <= sti <= length(wing.strips))
        @assert(n >= 0)
        set_wing_to_single_bv_fn!(fil_wing, wing, fil_locs, n, sti)
        ind_vel_n = x->ind_vel_external(x) + ind_vel(fil_wing, x)    # Is this correct?
        vf = get_particle_vorticity_function(
            particles, k_sind, dt, ind_vel_n)
        particle_vorts = vf(wing_to_bv_vector(fil_wing), old_bound_vorticities)
        for i = 1 : length(particle_vorts)
            particles[i].vorticity =
                unit(particles[i].vorticity) * particle_vorts[i]
        end
        # The purely 2D bit on the wing.
        u = (dx, dy)-> (dx<=eps(Float64) && dy<=eps(Float64) ? 0 : dy / (2 * pi * (dy^2 + dx^2)))
        v = (dx, dy)-> (dx<=eps(Float64) && dy<=eps(Float64) ? 0 :- dx / (2 * pi * (dy^2 + dx^2)))
        svort_locs = surf_fn.(sti, fil_locs)
        density_func = vorticity_density_x_fn(n)
        vort = lump_vorticities(wing.strips[sti], density_func, fil_locs)
        function ind_vel_2d(x :: ThreeDVector)
            dx = map(y->x.x - y.x, svort_locs)
            dy = map(y->x.z - y.z, svort_locs)
            velu = mapreduce(x->u(x[1], x[2]) * x[3], +, 0.0, zip(dx, dy, vort))
            velv = mapreduce(x->v(x[1], x[2]) * x[3], +, 0.0, zip(dx, dy, vort))
            return ThreeDVector(0, velu, velv)
        end
        ind_vel_pw = x-> ind_vel(particles, x) + ind_vel(fil_wing, x) - ind_vel_2d(x)
        return ind_vel_pw
    end

    nex, ney, nez = nocamber_normal_splines(wing)
    normal_dir_fn = x->ThreeDVector(nex(x), ney(x), nez(x))
    deta_dx_dot_c_fn = get_surface_detadx_dot_c_fn(wing)
    surf_fn = get_surface_fn(wing)
    mtrx = downwash_2D_fourier_integrals_all_strips(
        normal_dir_fn, deta_dx_dot_c_fn, surf_fn,
        length(wing.strips), n_max, ref_vel, reset, fil_locs )
    return mtrx
end


function compute_new_particles_fourier_integral_matrix(
    wing :: StripDefinedWing,
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    fil_locs :: Vector{Float64},
    old_bound_vorticities :: Vector{Float64},
    ref_vel :: Float64,
    new_particles :: ThreeDVortexParticleSet,
    kelvin_particle_indexes :: Vector{Int64},
    dt :: Float64,
    ind_vel_external :: Function,
    n_max :: Int64
    )
    @assert(length(wing.strips) > 0)
    @assert(length(wing.strips) == length(fil_wing.filaments_yp))
    @assert(isfinite(ref_vel))
    @assert(length(new_particles) ==
                length(kelvin_particle_indexes) + length(wing.strips) * 2 + 1)
    @assert(all(-1 .< fil_locs .< 1))
    @assert(n_max > 0)
    @assert(dt > 0)

    ind_vel_np :: Function = x->x # Induced vel excluding new particles. [placeholder for scope]
    ind_vel_pw :: Function = x->ind_vel(fil_wing, x) + ind_vel(particles, x) # Induced vel due to new particles and wing. [placeholder for scope]
    function reset(sti :: Int64, n :: Int64)
        @assert(1 <= sti <= length(wing.strips))
        @assert(n >= 0)
        set_wing_to_single_bv_fn!(fil_wing, wing, fil_locs, n, sti)
        ind_vel_n = x->ind_vel_external(x) + ind_vel(fil_wing, x)    # Is this correct?
        vf = get_particle_vorticity_function(
            particles, k_sind, dt, ind_vel_n)
        particle_vorts = vf(wing_to_bv_vector(fil_wing), old_bound_vorticities)
        for i = 1 : length(particle_vorts)
            particles[i].vorticity =
                unit(particles[i].vorticity) * particle_vorts[i]
        end
        ind_vel_pw = x-> ind_vel(particles, x)
        return ind_vel_pw
    end
    nex, ney, nez = nocamber_normal_splines(wing)
    normal_dir_fn = x->ThreeDVector(nex(x), ney(x), nez(x))
    deta_dx_dot_c_fn = get_surface_detadx_dot_c_fn(wing)
    surf_fn = get_surface_fn(wing)
    mtrx = downwash_2D_fourier_integrals_all_strips(
        normal_dir_fn, deta_dx_dot_c_fn, surf_fn,
        length(wing.strips), n_max, ref_vel, reset, fil_locs )
    return mtrx
end


function solve_new_vortex_particle_vorticities_and_assign!(
    untransformed_wing :: StripDefinedWing,
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    kinem :: ThreeDCoordinateTransform,
    fil_locs :: Vector{Float64},
    old_bound_vorticities :: Vector{Float64},
    ref_vel :: Float64,
    new_particles :: ThreeDVortexParticleSet,
    kelvin_particle_indexes :: Vector{Int64},
    dt :: Float64,
    ind_vel_external :: Function,
    n_fourier_terms :: Int64
    )
    @assert(all(-1 .< fil_locs .< 1))
    wing = transform_ThreeDVectors(x->kinem(x), untransformed_wing)
    nex, ney, nez = nocamber_normal_splines(wing)
    normal_dir_fn = x->ThreeDVector(nex(x), ney(x), nez(x))
    deta_fn = get_surface_detadx_dot_c_fn(wing)
    surf_fn = get_surface_fn(wing)

    # Downwash matrix due to the new wing / vortex particle combo.
    p_w_mtrx = compute_new_particles_and_wing_fourier_integral_matrix(
        wing, fil_wing, fil_locs, old_bound_vorticities,
        ref_vel, new_particles, kelvin_particle_indexes,
        dt, ind_vel_external, n_fourier_terms)

    # The downwash due to the wing's velocity.
    w_vel_vec = wing_velocity_to_fourier_integrals_vector(
        normal_dir_fn, deta_fn, untransformed_wing, kinem,
        length(wing.strips), n_fourier_terms,
        ref_vel, fil_locs )

    # The downwash due to the wing on itself in a 2D sense
    #=dw_2d = downwash_2D_pure_2D_self_fourier_integrals_matrix(
            wing, filament_positions, n_fourier_terms, ref_vel)=#

    # The downwash due to the wake and free stream
    dw_wake = downwash_2D_fourier_integrals_all_strips(
        normal_dir_fn, deta_fn, surf_fn, ind_vel_external, length(wing.strips),
        n_fourier_terms, ref_vel, fil_locs)

    fourier_terms = (p_w_mtrx - eye(size(p_w_mtrx)[1])) \ -(w_vel_vec + dw_wake)
    strip_ft = Vector{Vector{Float64}}(length(wing.strips))
    for s = 1 : length(wing.strips)
        idxs = fourier_term_linear_index.(
            n_fourier_terms, length(wing.strips), 0:n_fourier_terms, s)
        strip_ft[s] = fourier_terms[idxs]
    end
    set_wing_to_fourier_set!(fil_wing, wing, filament_positions, strip_ft)
    bvs = wing_to_bv_vector(fil_wing)
    ind_vel_fn = x->ind_vel_external(x) + ind_vel(fil_wing, x)
    bvf = get_particle_vorticity_function(new_particles,
        kelvin_particle_indexes, dt, ind_vel_fn)
    particle_vorts = bvf(bvs, old_bvs)
    for i = 1 : length(particles)
        particles[i].vorticity = particle_vorts[i] * unit(particles[i].vorticity)
    end
    return strip_ft
end


function pressure_distribution(
    current_fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    old_fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    current_wing :: StripDefinedWing,
    filament_pos :: Vector{Float64},
    dt :: Float64,
    ind_vel_fn :: Function,
    fourier_coeffs :: Vector{Vector{Float64}},
    density :: Float64
    )
    # We evaluate pressure at the centre of each vortex filament.
    ns = length(current_wing.strips)
    spos = Vector{Float64}(0.75 : 0.5 : ns + 0.25)  # Places on the span.
    cpos = filament_pos                             # Places on the chord
    cidx = Vector{Int64}(length(spos) * length(cpos)) # The idx on corresponding to cpos
    sidx = Vector{Int64}(length(spos) * length(cpos)) # The idx on corresponding to spos
    locations = Vector{ThreeDVector}(length(spos) * length(cpos)) # Pressure eval locs
    velocities = Vector{ThreeDVector}(length(spos) * length(cpos))  # Flow velocity at locs
    spanwise_tang = Vector{ThreeDVector}(length(spos) * length(cpos))   # tangent in s dir
    chordwise_tang = Vector{ThreeDVector}(length(spos) * length(cpos))  # tangent in c dir
    param_loc_s = Vector{Float64}(length(spos) * length(cpos))  # parametric spanwise loc
    param_loc_c = Vector{Float64}(length(spos) * length(cpos))  # parametric chordwise loc
    pressure = Vector{Float64}(length(spos) * length(cpos)) # Pressure at locations
    dvortdc = Vector{Float64}(length(spos) * length(cpos))  # Vorticity deriv chordwise
    dvortds = Vector{Float64}(length(spos) * length(cpos))  # Vorticity deriv spanwise

    for i = 1 : length(spos)
        for j = 1 : length(cpos)
            idx = (i - 1) * length(cpos) + j
            param_loc_s[idx] = spos[i]
            param_loc_c[idx] = cpos[j]
            sidx[idx] = i
            cidx[idx] = j
        end
    end

    slant_cors = Vector{Vector{Float64}}(length(spos))
    for i = 1 : ns
        scym, scyp = slant_correction_factors(
            current_wing, current_fil_wing, cpos, i)
        slant_cors[i * 2 - 1] = scym
        slant_cors[i * 2] = scyp
    end

    surf_fn = get_surface_fn(current_wing)
    for i = 1 : (length(spos) * length(cpos))
        sp = param_loc_s[i] # Span parametric coord
        stripidx = Int64(round(sp))
        cp = param_loc_c[i] # Chord parametric coord
        locations[i] = surf_fn(sp, cp)
        # Annoyingly the simple thing to do here is numerical differentiation.
        spanwise_tang[i] = (surf_fn(sp+0.0001, cp) -surf_fn(sp-0.0001, cp)) / 0.02
        chordwise_tang[i] = (surf_fn(sp, cp+0.0001) -surf_fn(sp, cp-0.0001)) / 0.02
        velocities[i] = ind_vel_fn(locations[i])
        dvortdc = vorticity_density_x_derivative(fourier_coeffs[stripidx], cp) *
            slant_cors[stripidx][cidx[i]]
        # Averaging of the chordwise vorticity (remembering dirs of filament
        # is always TE to LE.)
        # Also the rate of change of strenght of the vortex ring wrt time.
        fcidx = cidx[i] == length(cpos) ? cidx[i] - 1 : cidx[i]
        dvortds = (
            fil_wing.filaments_cw[sidx[i]][fcidx].vorticity
            + fil_wing.filaments_cw[sidx[i] + 1][fcidx].vorticity
            ) / 2
        if sp < round(sp)
            dvdt = mapreduce(x->x[1].vorticity - x[2].vorticity, +, 0.0,
                            zip(current_fil_wing.filaments_ym[stripidx],
                                old_fil_wing.filaments_ym[stripidx])
                            ) / dt
        else
            dvdt = mapreduce(x->x[1].vorticity - x[2].vorticity, +, 0.0,
                            zip(current_fil_wing.filaments_yp[stripidx],
                                old_fil_wing.filaments_yp[stripidx])
                            ) / dt
        end
        pressure[i] = density *
            ( dvdt)
        #=pressure[i] = density *
            (dot(velocities[i], spanwise_tang[i] * dvortds
                + chordwise_tang[i] * dvortdc) + dvdt)=#
    end
    println(pressure)
    pressure_field = Spline2D(param_loc_c, param_loc_s,
                                    pressure, kx=1, ky=1, s=length(param_loc_c))
    return pressure_field
end

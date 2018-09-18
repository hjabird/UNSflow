#===============================================================================
calcs.jl

Calculations for lowOrder3D methods.
===============================================================================#

#===============================================================================
    ThreeDSurfSimple methods
------------------------------------------------------------------------------=#
function calc_a0a13d(surf::ThreeDSurfSimple)

    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)


    for i = 1:surf.nspan
        integ = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        -simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)

        for n = 1:surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + (nn*pi/(2*surf.AR))
            *(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ/pi))

        end
        rhs[i] = pi*sin(surf.psi[i])*surf.bc[i]/(2*surf.AR)
    end

    surf.bcoeff[:] = lhs \ rhs

    for i = 1:surf.nspan
        surf.a03d[i] = 0
        surf.aterm3d[1,i] = 0
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)

        for n = 1:surf.nspan
            nn = 2*n - 1
            surf.a03d[i] -= real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])*(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ0/pi)

            surf.aterm3d[1,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ1


        end
    end
end

function calc_a2toan3d(surf::ThreeDSurfSimple)
    for ia = 2:surf.naterm
        for i = 1:surf.nspan
            surf.aterm3d[ia,i] = 0
            integ = simpleTrapz(surf.s2d[i].cam_slope.*cos.(ia*surf.s2d[i].theta), surf.s2d[i].theta)

            for n = 1:surf.nspan
                nn = 2*n - 1
                surf.aterm3d[ia,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ
            end
        end
    end
end
 # END ThreeDSurfSimple methods ===============================================#

#===============================================================================
      ThreeDVortexParticleSet functions

      Initial code: HJAB 2018
------------------------------------------------------------------------------=#
"""
Define the vortex particle interaction kernal.

Assign a kernal such as threed_singularity_kernels() to the vortex particle set.
"""
function set_particle_interaction_function!(
    particle_set :: ThreeDVortexParticleSet,
    threed_kernal :: Function
    )

    g, f = threed_kernal()
    particle_set._reduction_factor_fn = g
    particle_set._vorticity_fraction_fn = f
    return
end

function mutual_ind(
    particle_set :: ThreeDVortexParticleSet
    )

    return mutual_ind(particle_set.particles, particle_set._reduction_factor_fn,
        _vorticity_fraction_fn)
end

function ind_vel(
	particle_set :: ThreeDVortexParticleSet,
	mes_pos :: ThreeDVector
	)
	return ind_vel(particle_set.particles, mes_pos,
		particle_set._reduction_factor_fn)
end


""" Sums any externally induced velocity or rate of change of vorticity
to the existing.

external_velocity_influence is a function that is takes a single
ThreeDVortexParticle argument representing a cartesian coordinate and returns a
ThreeDVector representing a an induced velocity.

external_rate_of_change_of_vorticity_influence is a function that is takes a
single ThreeDVortexParticle argument representing a cartesian coordinate and
returns a ThreeDVector representing the induced rate of change of vorticity.
"""
function external_ind!(
    particle_set :: ThreeDVortexParticleSet,
    external_velocity_influence :: Function,
    external_rate_of_change_of_vorticity_influence :: Function
    )

    n = size(particle_set.particles)[1]
    for i = 1 : n
        particle_set.particles[i].velocity +=
            external_velocity_influence(particle_set.particles[i].coord)
        particle_set.particles[i].vorticity_time_derivative +=
            external_rate_of_change_of_vorticity_influence(particle_set.particles[i])
    end
    return particle_set
end

function integration_step!(
    particle_set :: ThreeDVortexParticleSet,
    dt :: Float64
    )
    integration_step!(particle_set.particles, dt)
    return particle_set
end

"""
Compute velocities and rate of change of vorticities, and convect vortex
particles subject to external influence according to forward Euler method.

external_velocity_influence is a function that is takes a single
ThreeDVortexParticle argument representing a cartesian coordinate and returns a
ThreeDVector representing a an induced velocity.

external_rate_of_change_of_vorticity_influence is a function that is takes a
single ThreeDVortexParticle argument representing a cartesian coordinate and
returns a ThreeDVector representing the induced rate of change of vorticity.
"""
function euler_forward_step!(
    particle_set :: ThreeDVortexParticleSet,
    external_velocity_influence :: Function,
    external_rate_of_change_of_vorticity_influence :: Function,
    dt :: Float64
    )

    if length(particle_set) == 0
        return particle_set
    end
    mutual_ind(particle_set.particles, particle_set._reduction_factor_fn,
        particle_set._vorticity_fraction_fn)
    external_ind!(particle_set, external_velocity_influence,
        external_rate_of_change_of_vorticity_influence)
    integration_step!(particle_set.particles, dt)
    return particle_set
end

"""
Compute velocities and rate of change of vorticities, and convect vortex
particles subject to external influence according to explicit midpoint method.

external_velocity_influence is a function that is takes a single
ThreeDVortexParticle argument representing a cartesian coordinate and returns a
ThreeDVector representing a an induced velocity.

external_rate_of_change_of_vorticity_influence is a function that is takes a
single ThreeDVortexParticle and
returns a ThreeDVector representing the induced rate of change of vorticity.
"""
function explicit_midpoint_step!(
    particle_set :: ThreeDVortexParticleSet,
    external_velocity_influence :: Function,
    external_rate_of_change_of_vorticity_influence :: Function,
    dt :: Float64
    )


    mutual_ind(
        particle_set.particles,
        particle_set._reduction_factor_fn,
        particle_set._vorticity_fraction_fn
        )
    external_ind!(particle_set, external_velocity_influence,
        external_rate_of_change_of_vorticity_influence)
    particle_set_cpy = deepcopy(particle_set)
    integration_step!(particle_set_cpy, dt / 2.)

    mutual_ind(particle_set.particles,
        particle_set._reduction_factor_fn,
        particle_set._vorticity_fraction_fn
        )
    external_ind!(particle_set, external_velocity_influence,
        external_rate_of_change_of_vorticity_influence)
    for (particle, particlecpy) in zip(particle_set.particles,
            particle_set_cpy.particles)
        particle.velocity = particlecpy.velocity
        particle.vorticity_time_derivative =
            particlecpy.vorticity_time_derivative
    end
    integration_step!(particle_set, dt)
    return particle_set
end
#= END ThreeDVortexParticleSet functions ======================================#

#===============================================================================
      ThreeDStraightVortexFilament functions

      Initial code: HJAB 2018
------------------------------------------------------------------------------=#
function ind_vel(
    inducing_filament :: ThreeDStraightVortexFilament,
    measurement_loc :: ThreeDVector
    )
    r0 = inducing_filament.end_coord - inducing_filament.start_coord
    r1 = measurement_loc - inducing_filament.start_coord
    r2 = measurement_loc - inducing_filament.end_coord
    if(abs(r1) <= eps(Float64) || abs(r2) <= eps(Float64))
        return ThreeDVector(0,0,0)
    end
    # From Katz & Plotkin, Eq(2.72), pg41
    term1 = inducing_filament.vorticity / (4 * pi)
    term2n = cross(r1, r2)
    term2d = abs(cross(r1, r2)) ^ 2
    term3 = dot(r0, unit(r1) - unit(r2))
    vel =  term1 * (term2n / term2d) * term3
    vel = (!isfinite(vel) ? ThreeDVector(0,0,0) : vel)
    return vel
end

function ind_vel(
    inducing_filament :: ThreeDStraightVortexFilament,
    measurement_locs :: Vector{ThreeDVector}
    )
    return map(x->ind_vel(inducing_filament, x), measurement_locs)
end

function ind_vel(
    inducing_filaments :: Vector{ThreeDStraightVortexFilament},
    measurement_loc :: ThreeDVector
    )
    return map(x->ind_vel(x, measurement_loc), inducing_filaments)
end

function ind_dvortdt(
    induced_particle :: ThreeDVortexParticle,
    inducing_filament :: ThreeDStraightVortexFilament
    )
    r0 = inducing_filament.end_coord - inducing_filament.start_coord
    r1 = induced_particle.coord - inducing_filament.start_coord
    r2 = induced_particle.coord - inducing_filament.end_coord
    # Notes, HJAB, Book 4, pg.42-pg.43
    term1 = inducing_filament.vorticity / ( 4 * pi)
    term211 = cross(induced_particle.vorticity, r0) / (abs(cross(r1, r0))^2)
    term2121 = dot(r0, r1) / abs(r1)
    term2122 = -dot(r0, r2) / abs(r2)
    term221 = 3.0 * induced_particle.vorticity / abs(r0)
    term2221 = abs(cross(r0, r1)) / abs(r1)
    term2222 = -abs(cross(r0, r2)) / abs(r2)
    return term1 * (term211 * (term2121 + term2122) +
        term221 * (term2221 + term2222))
end

function ind_dvortdt(
    induced_particles :: Vector{ThreeDVortexParticle},
    inducing_filament :: ThreeDStraightVortexFilament
    )

    return map(x->ind_dvortdt(x, inducing_filament), induced_particle)
end

function ind_dvortdt(
    induced_particles :: ThreeDVortexParticle,
    inducing_filament :: Vector{ThreeDStraightVortexFilament}
    )

    return map(x->ind_dvortdt(induced_particles, x), inducing_filament)
end
#= END ThreeDStraightVortexFilament functions =================================#

#===============================================================================
      ThreeDVortexRing functions

      Initial code: HJAB 2018
------------------------------------------------------------------------------=#
function ind_vel(
    inducing_ring :: ThreeDVortexRing,
    measurement_loc :: ThreeDVector
    )
    return ind_vel(Vector{ThreeDStraightVortexFilament}(inducing_ring),
                                                            measurement_loc)
end

function ind_dvortdt(
    induced_particle  :: ThreeDVortexParticle,
    inducing_ring :: ThreeDVortexRing
    )
    return ind_dvortdt(Vector{ThreeDStraightVortexFilament}(inducing_ring),
                                                            measurement_loc)
end

function area(vortex_ring :: ThreeDVortexRing)
    a = vortex_ring
    # Gauss-Legendre quadrature
    x = 1 / sqrt(3)
    xs = [x, -x, x, -x]
    ys = [x, x, -x, -x]
    ws = [1, 1, 1, 1]
    # Derivatives of a bilinear surface.
    dfdx = (x,y)-> 0.25 * ((y-1)*a.c1 - (y-1)*a.c2 + (y+1)*a.c3 - (y+1)*a.c4)
    dfdy = (x,y)-> 0.25 * ((x-1)*a.c1 - (x+1)*a.c2 + (x+1)*a.c3 - (x-1)*a.c4)
    integrand = (x, y)->abs(cross(dfdx(x, y), dfdy(x,y)))
    integral = mapreduce(x->x[1] * integrand(x[2],x[3]),
                        +, 0.0, zip(ws, xs, ys))
    return integral
end

function normal(vortex_ring :: ThreeDVortexRing)
    return cross(tangent_dir1(vortex_ring), tangent_dir2(vortex_ring))
end

function tangent_dir1(vortex_ring :: ThreeDVortexRing)
    a = vortex_ring
    # Based on bilinear interpolation at (0,0)
    return 0.25 * (-a.c1 + a.c2 + a.c3 - a.c4)
end

function tangent_dir2(vortex_ring :: ThreeDVortexRing)
    a = vortex_ring
    # Based on bilinear interpolation at (0,0)
    return 0.25 * (-a.c1 - a.c2 + a.c3 + a.c4)
end

function centre(vortex_ring :: ThreeDVortexRing)
    a = ThreeDVortexRing
    return (a.c1 + a.c2 + a.c3 + a.c4) / 4
end

function pressure_difference(
    vortex_ring :: ThreeDVortexRing,
    ind_vel_fn :: Function,
    density :: Float64,
    old_circ_strength :: Float64,
    dt :: Float64,
    circ_ip :: Float64, # Adjacent ring circulation strengths in i / j plus dir.
    circ_jp :: Float64
    )
    @assert(density >= 0)
    @assert(isfinite(old_circ_strength))
    @assert(isfinite(circ_ip))
    @assert(isfinite(circ_jp))
    @assert(isfinite(vortex_ring.c1))
    # Based on Katz & Plotkin 13.12 F (pg. 427)
    a = vortex_ring
    x = centre(vortex_ring)
    tangent1 = tangent_dir1(vortex_ring)
    tangent2 = tangent_dir2(vortex_ring)
    c_len = abs(a.c2 - a.c1 + a.c3 - a.c4) / 2
    b_len = abs(a.c4 + a.c3 - a.c2 - a.c1) / 2
    # B differences
    di = (circ_ip - a.strength) / c_len
    dj = (circ_jp - a.strength) / b_len
    vel = ind_vel_fn(x)
    p = density * (dot(vel, tagent1 * di + tangent2 * dj)) +
        (a.strength - old_circ_strength) / dt
    return p
end

function force(
    vortex_ring :: ThreeDVortexRing,
    ind_vel_fn :: Function,
    density :: Float64,
    old_circ_strength :: Float64,
    dt :: Float64,
    circ_ip :: Float64, # Adjacent ring circulation strengths in i / j plus dir.
    circ_jp :: Float64)
    dp = pressure_difference(vortex_ring, ind_vel_fn, density,
        old_circ_strength, dt, circ_ip, circ_jp)
    f = dp * area(vortex_ring)
    return f
end
#= END ThreeDVortexRing functions =============================================#


#===============================================================================
      ThreeDSpanwiseFilamentWingRepresentation functions

      Initial code: HJAB 2018
------------------------------------------------------------------------------=#
function transform_wing_location(
    wing :: ThreeDSpanwiseFilamentWingRepresentation,
    kinematics :: KinemDef3D
    )
    for ic = 1 : wing.n_chords
        for ifil = 1 : n_filaments_per_chord[ic]
            mod!(nwing.filaments_ym[i][j].start_coord)
            mod!(nwing.filaments_ym[i][j].end_coord)
            mod!(nwing.filaments_yp[i][j].start_coord)
            mod!(nwing.filaments_yp[i][j].end_coord)
            mod!(nwing.filaments_cw[i][j].start_coord)
            mod!(nwing.filaments_cw[i][j].end_coord)
        end
    end
    return nwing
end

"""
Generate corrections for the slant of vortex filaments in comparason to the
plane in which the vorticities are being computed.
"""
function slant_correction_factors(
    wing :: StripDefinedWing,
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    filament_locs :: Vector{Float64}, # in [-1, 1]
    strip_idx :: Int64
    )
    @assert(strip_idx > 0)
    @assert(strip_idx <= length(fil_wing.filaments_yp))
	@assert(length(fil_wing.filaments_yp) == length(fil_wing.filaments_ym))
    @assert(all(abs.(filament_locs) .<= 1.0))

    c = chord(wing.strips[strip_idx])
    nsx, nsy, nsz = nocamber_normal_splines(wing)
    normal = ThreeDVector(map(f->f(Float64(strip_idx)), [nsx, nsy, nsz]))
    spanwise = unit(cross(c, normal))

    corrections_ym = zeros(size(filament_locs)[1])
    corrections_yp = zeros(size(filament_locs)[1])
    function corr(d :: ThreeDStraightVortexFilament)
		a = d.start_coord
		b = d.end_coord
        return 1. / dot(unit(b-a), spanwise)
    end
    for i = 1 : size(filament_locs)[1]
        corrections_ym[i] = abs.(corr(fil_wing.filaments_ym[strip_idx][i]))
        corrections_yp[i] = abs.(corr(fil_wing.filaments_yp[strip_idx][i]))
    end
    @assert(all(corrections_ym .>= 1.0 * (1 - eps(Float64))))
    @assert(all(corrections_yp .>= 1.0 * (1 - eps(Float64))))
    return corrections_ym, corrections_yp
end

"""
Helmholtz filaments (the ones that are due to vortex filments being never
ending, present even in steady state) are shed at all 2 * ns + 1 junctions
in the filament wing. This returns the strength of the filament
"""
function helmholtz_vortex_strenghts(
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation
    )
    ns = fil_wing.n_strips
    strengths = zeros(ns * 2 + 1)
    f_yp = zeros(ns)
    f_ym = zeros(ns)
    # On the filament_yp <-> filament_ym interface / tips
    for i = 1 : ns
        f_yp = mapreduce(x->x.vorticity, +, 0.0, fil_wing.filaments_yp[i])
        f_ym = mapreduce(x->x.vorticity, +, 0.0, fil_wing.filaments_ym[i])
    end
    strengths[1] = f_ym[1]
    for i = 1 : ns - 1
        strengths[i * 2 + 1] = f_ym[i + 1] - f_yp[i]
    end
    strengths[end] = f_yp[end]
    # On the filament_ym <-> filament_yp interface (the defined chords)
    for i = 1 : ns
        strengths[i * 2] = f_yp[i] - f_ym[i]
    end
    return strengths
end

""" Helmholtz filaments (the ones that are due to vortex filments being never
ending, present even in steady state) are shed at all 2 * ns + 1 junctions
in the filament wing. This returns the strength of the filament DUE TO
ONLY A SINGLE STRIP'S vorticity.
"""
function helmholtz_vortex_strenghts(
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    strip_idx :: Int64
    )
    strengths = zeros(3)
    f_yp = mapreduce(x->x.vorticity, +, 0.0, fil_wing.filaments_yp[strip_idx])
    f_ym = mapreduce(x->x.vorticity, +, 0.0, fil_wing.filaments_ym[strip_idx])
    strengths[1] = f_ym
    strengths[2] = f_yp - f_ym
    strengths[3] = -f_yp
    return strengths
end


function new_vortex_particle_locations(
    wing :: StripDefinedWing,
    old_vortex_particles_positions :: Vector{ThreeDVector},
    span_shed_loc :: Vector{Float64}
    )
    @assert(all(span_shed_loc .>= 0.5))
    @assert(all(span_shed_loc .<= length(wing.strips) + .5))
    @assert(size(span_shed_loc)
        == size(old_vortex_particles_positions))
    oppos = old_vortex_particles_positions
    tex, tey, tez = te_spline(wing)
    pos = span_shed_loc
    tepos = map(x->ThreeDVector(tex(x), tey(x), tez(x)), pos)
    vpos = map(x->tepos[x] + (oppos[x] - tepos[x]) / 3.0, 1 : size(pos)[1])
    return vpos
end

function helmholtz_vortex_particle_locations(
    wing :: StripDefinedWing,
    old_helmholtz_vortex_particles_positions :: Vector{ThreeDVector}
    )
	@assert(length(old_helmholtz_vortex_particles_positions) ==
		length(wing.strips) * 2 + 1)
    pos = Vector{Float64}(.5 : .5 : Float64(length(fil_wing.filaments_ym)) + .5)
    return new_vortex_particle_locations(wing,
        old_helmholtz_vortex_particles_positions, pos)
end

function initial_helmholtz_particles_positions(
    wing :: StripDefinedWing,
    free_stream :: Function,
    dt :: Float64
    )

    if isfinite(dt) != true
        error("Expected finite value for time step.")
    end
    tex, tey, tez = te_spline(wing)
    const n_strips = length(wing.strips)
    const n_p = 2 * n_strips + 1
    initial_locs = Vector{ThreeDVector}(n_p)
    function loc(x :: Float64)
        locn = ThreeDVector(tex(x), tey(x), tez(x))
        dx = free_stream(locn) * dt * 0.5
        if isfinite(dx) != true
            error("Evaluated free_stream_vel * dt as non-finite!")
        end
        return locn + dx
    end
	for i = 0.5 : 0.5 : n_strips + 0.5
		initial_locs[Int64(i * 2)] = loc(i)
	end
    return initial_locs
end

"""
Find the locations on the TE to shed kelvin vortex particles.

This function returns a vector of floats that relate to the TE_spline evaluation
positions. This allows relation to the different sections on the TE.
"""
function kelvin_particles_span_shedding_locations(
    wing :: StripDefinedWing,
    max_particle_size :: Float64
    )

    ns = length(wing.strips)
    tex, tey, tez = te_spline(wing)
    f = x->ThreeDVector(tex(x), tey(x), tez(x))
    function subsection(strip :: Int64)
		# Places to evaluate the tailing edge spline:
        tm = Vector{Float64}(strip - 0.5 : 0.05: strip)
        tp = Vector{Float64}(strip : 0.05: strip + 0.5)
		# Eval spline
        ftm = f.(tm)
        ftp = f.(tp)
        # Compute the lengths between the adjacent points and sum.
        ftml = sum(abs.(ftm[2 : end] - ftm[1 : end - 1]))
        ftpl = sum(abs.(ftp[2 : end] - ftp[1 : end - 1]))
        # The number of particles we'll need
        n_ym = ceil(ftml / max_particle_size)
        n_yp = ceil(ftpl / max_particle_size)
        # The where to evaluate our spline to obtain te locations
        pmt = Vector{Float64}(linspace(tm[1], tm[end], n_ym + 1))
        pmt = 0.5 * (pmt[1:end-1] + pmt[2:end])
        ppt = Vector{Float64}(linspace(tp[1], tp[end], n_yp + 1))
        ppt = 0.5 * (ppt[1:end-1] + ppt[2:end])
		positions = vcat(pmt, ppt)
        return positions
    end

	positions = Vector{Float64}([])
    for i = 1 : ns
        positions = vcat(positions, subsection(i))
    end
    return positions
end

"""
Given the shedding locations of the kelvin particles, find the positions
for the first time step
"""
function initial_kelvin_particle_locations(
    shedding_locations :: Vector{ThreeDVector},
    free_stream :: Function,
    dt :: Float64
    )
    if isfinite(dt) != true
        error("Expected finite value for time step.")
    end
    new_locs = deepcopy(shedding_locations)
    for i = 1 : length(new_locs)
		dx = free_stream(new_locs[i]) * dt * 0.5
        new_locs[i] += free_stream(new_locs[i]) * dt * 0.5
        if isfinite(i) != true
            error("Expected finite valued locations")
        end
    end
    return new_locs
end

"""
Shed initial particles.

Particle have no vorticity.
"""
function shed_initial_particles(
	wing :: StripDefinedWing,
	free_stream :: Function,
	k_particle_shedding_locations :: Vector{Float64},
	dt :: Float64,
	particle_size :: Float64
	)
	tex, tey,tez = te_spline(wing)
	fte = x -> ThreeDVector(tex(x), tey(x), tez(x))
	k_locs = initial_kelvin_particle_locations(
		map(fte, k_particle_shedding_locations), free_stream, dt)
	h_locs = initial_helmholtz_particles_positions(wing, free_stream, dt)
	nk = length(k_locs)
	nh = length(h_locs)
	particles = Vector{ThreeDVortexParticle}(nk + nh)
	for i = 1 : nh
		particles[i] = ThreeDVortexParticle(
			h_locs[i],
			ThreeDVector(0,0,0),
			particle_size,
			ThreeDVector(0,0,0),
			ThreeDVector(0,0,0)
		)
	end
	for i = nh + 1 : nh + nk
		particles[i] = ThreeDVortexParticle(
			k_locs[i - nh],
			ThreeDVector(0,0,0),
			particle_size,
			ThreeDVector(0,0,0),
			ThreeDVector(0,0,0)
		)
	end
	return ThreeDVortexParticleSet(particles)
end

"""
Shed TE particles.
Returns particle vector. No vorticities assigned.
"""
function shed_particles(
	wing :: StripDefinedWing,
	k_particle_shedding_locations :: Vector{Float64},
	particle_size :: Float64,
	old_particles :: ThreeDVortexParticleSet
	)
	nh = length(wing.strips) * 2 + 1
	nk = length(old_particles) - nh
	particles = Vector{ThreeDVortexParticle}(nh + nk)
	plocs = Vector{ThreeDVector}(nh + nk)
	plocs[1:nh] =
		helmholtz_vortex_particle_locations(wing,
			map(x->x.coord, old_particles[1:nh]))
	plocs[nh + 1 : nh + nk] =
		new_vortex_particle_locations(wing,
			map(x->x.coord, old_particles[nh + 1: end]),
			k_particle_shedding_locations)
	for i = 1 : length(particles)
		particles[i] = ThreeDVortexParticle(
			plocs[i],
			ThreeDVector(0,0,0),
			particle_size,
			ThreeDVector(0,0,0),
			ThreeDVector(0,0,0)
			)
	end
	return ThreeDVortexParticleSet(particles)
end

"""
Convert the Vector{Float64} of locations to shed k_vortices to indexes
associating the vorticies to a bound vorticity section on the wing.
"""
function k_particle_shedding_locs_to_bv_index(
	s_loc :: Vector{Float64},
	n_strips :: Int64
	)
	@assert(all(s_loc .>= 0.5))
	@assert(all(s_loc .<= n_strips + 0.5))
	# s_loc goes from 0 to num_strips + 1.
	# bit associated with a strip is ns +- 0.5. -0.5 is y_minus side,
	# +0.5 is y_plus side. Except at the tips where it is -1 and +1...
	a = deepcopy(s_loc)
	b = Vector{Int64}(length(a))
	for i = 1 : length(a)
		c = round(a[i])
		b[i] = c > a[i] ? 2 * c - 1 : 2 * c
	end
	return b
end

"""
Obtain a bound vorticity vector from the wing.

An n strip wing will return a 2*n vector.
"""
function wing_to_bv_vector(
	fil_wing :: ThreeDSpanwiseFilamentWingRepresentation
	)
	@assert(length(fil_wing.filaments_yp)
									== length(fil_wing.filaments_ym))
	ns = length(fil_wing.filaments_ym)
	bvs = zeros(ns * 2)
	for i = 1 : ns
		bvs[2*i-1] = mapreduce(x->x.vorticity, +, 0.0, fil_wing.filaments_ym[i])
		bvs[2*i] = mapreduce(x->x.vorticity, +, 0.0, fil_wing.filaments_yp[i])
	end
	@assert(length(bvs) == ns * 2)
	@assert(all(isfinite.(bvs)))
	return bvs
end

"""
Returns a function to obtain the intensity of particle vorticities
from a bound vorticity vector
"""
function get_particle_vorticity_function(
	particles :: ThreeDVortexParticleSet, # ordered nh, nk
	k_particle_to_bv_index_mapping :: Vector{Int64},
	dt :: Float64,
	ind_vel :: Function
	)

	nk = length(k_particle_to_bv_index_mapping)
	nh = length(particles) - nk
	if all(k_particle_to_bv_index_mapping .< nh) != true
		error("Some index mappings are wrong - are the function inputs right?")
	end
	if any(k_particle_to_bv_index_mapping .== nh - 1) != true
		error(string("It looks like at least the final segment of bound",
		" vorticity has no kelvin particles? ", nh, " Helmholtz particles",
		" and ", nk, " kelvin particles with max bv idx of ",
		maximum(k_particle_to_bv_index_mapping), "."))
	end
	bv_mat = zeros(nk + nh, nh - 1)
	for i = 1 : nh
		loc = particles[i].coord
		vel = ind_vel(loc)
		particles[i].velocity = vel
		particles[i].vorticity = unit(vel)
		if i < nh
			bv_mat[i, i] = abs(vel) * dt
		end
		if i > 1
			bv_mat[i, i - 1] = -abs(vel) * dt
		end
	end
	# Get dx for the k particles
	dx = Vector{ThreeDVector}(nk)
	x = map(x->x.coord, particles[nh + 1: nk + nh])
	dx[2 : end - 1] = (x[3 : end] .- x[1 : end - 2]) / 2.
	dx[1] = x[2] - x[1]
	dx[end] = x[end] - x[end - 1]
	for i = nh + 1 : nh + nk
		loc = particles[i].coord
		particles[i].vorticity = unit(dx[i - nh])
		idx = k_particle_to_bv_index_mapping[i - nh]
		bv_mat[i, idx] = abs(dx[i - nh])
	end
	function bv_to_vort_mult(
		bvs :: Vector{Float64}, old_bvs :: Vector{Float64})
		@assert(length(bvs) == length(old_bvs) == size(bv_mat)[2])
		r = zeros(nk + nh)
		r[1:nh] = bv_mat[1:nh, :] * bvs
		r[nh + 1 : end] = bv_mat[nh + 1 : nh + nk, :] * (bvs - old_bvs)
		return r
	end
	return bv_to_vort_mult
end


function ind_vel(
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    measurement_loc :: ThreeDVector,
    strip_idx :: Int64
    )
    @assert(0 < strip_idx <= length(fil_wing.filaments_ym))
    return sum(ind_vel(fil_wing.filaments_ym[strip_idx], measurement_loc)) +
        sum(ind_vel(fil_wing.filaments_yp[strip_idx], measurement_loc)) +
        sum(ind_vel(fil_wing.filaments_cw[strip_idx * 2], measurement_loc)) +
        sum(ind_vel(fil_wing.filaments_cw[strip_idx * 2], measurement_loc)) +
        sum(ind_vel(fil_wing.filaments_cw[strip_idx * 2], measurement_loc))
end

function ind_vel(
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    measurement_loc :: ThreeDVector
    )
    strips = mapreduce(
        x->ind_vel(fil_wing, measurement_loc, x),
        +, ThreeDVector(0,0,0), 1:length(fil_wing.filaments_ym))
    chords = mapreduce(
        x->sum(ind_vel(x, measurement_loc)),
        +, ThreeDVector(0,0,0), vec(fil_wing.filaments_cw)
    )
    return strips #+ chords
end

function ind_dvortdt(
    induced_particle :: ThreeDVortexParticle,
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    strip_idx :: Int64
    )

    return sum(
        ind_dvortdt(induced_particle, fil_wing.filaments_ym[strip_idx]) +
        ind_dvortdt(induced_particle, fil_wing.filaments_yp[strip_idx]))
end

function ind_dvortdt(
    induced_particle :: ThreeDVortexParticle,
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation
    )


    fils = convect(Vector{ThreeDStraightVortexFilament}, fil_wing)
    effect = mapreduce(
        x->ind_dvortdt(induced_particle, x),
        +, ThreeDVector(0,0,0), fils)
    return effect
end
#= END ThreeDSpanwiseFilamentWingRepresentation functions =====================#

#===============================================================================
      VortexParticleWakeLUATATSolution functions

      Initial code: HJAB 2018
------------------------------------------------------------------------------=#
function bound_vorticity(sim :: VortexParticleWakeLAUTATSolution)
    n_strips = size(sim.wing.strips)[1]
    bound_vs = zeros(n_strips)
    for i = 1 : n_strips
        mult = dot(chord(sim.wing.strips[i]), sim.free_stream_velocity) * pi
        bound_vs[i] = mult * (sim.fourier_terms[1] + fourier_terms[2] / 2)
    end
    return bound_vs
end

function vorticity_density(
    sim :: VortexParticleWakeLAUTATSolution,
    theta_pos :: Float64,
    idx :: Int64)

    @assert(idx <= size(sim.wing.strips)[1])
    @assert(idx > 0)

    mult = dot(unit!(chord(sim.wing.strips[i])), sim.free_stream_velocity) * 2
    a0 = sim.fourier_terms[idx][1] * (1 + cos(theta_pos)) / sin(theta_pos)
    an = mapreduce(x->sim.fourier_terms[x + 1] * sin(theta_pos * (x)),
        +, 0.0, 1 : sim.n_fourier_terms - 1)
    vort_d = mult * (a0 + an)
    return vort_d
end

function vorticity_density(
    sim :: VortexParticleWakeLAUTATSolution,
    theta_pos :: Float64
    )
    n_strips = size(sim.wing.strips)[1]
    vort_d = zeros(n_strips)
    for i = 1 : n_strips
        vort_d = vorticity_density(sim, theta_pos, i)
    end
    return vort_d
end

"""
Builds a ThreeDSpanwiseFilamentWingRepresentation from geometry defined by
StripDefinedWing.

chord_positions the locations that define the positioning of filaments in
[-1,1].

No vorticity is defined or kinematics applied.
"""
function build_vortex_filament_wing_geometry(
    wing :: StripDefinedWing,
    chord_positions :: Vector{Float64}  # In [-1, 1]
    )
    if all(abs.(chord_positions) .<= 1.0) != true
        error(string("Chord positions were not in [-1, 1] where -1 is leading",
            "edge and 1 is trailing edge."))
    end
    n_strips = length(wing.strips)
    if n_strips <= 0
        error("Wing must have a positive (nonzero) number of strips.")
    end
    n_fil = length(chord_positions)
    @assert(n_fil > 1)
    fil_wing = ThreeDSpanwiseFilamentWingRepresentation(n_strips, n_fil)
    surf_func = get_surface_fn(wing)
    for i in 1 : n_fil  # tip y_minus
        coord = surf_func(0.0, chord_positions[i])
        fil_wing.filaments_ym[1][i].start_coord = coord
        if i == 1
            fil_wing.filaments_cw[1][i].start_coord = coord
        elseif i < n_fil
            fil_wing.filaments_cw[1][i].start_coord = coord
            fil_wing.filaments_cw[1][i - 1].end_coord = coord
        else
            fil_wing.filaments_cw[1][i - 1].end_coord = coord
        end
    end
    for i_f in 1 : n_fil # In between the strips
        for i_s in 2 : n_strips
            coord = surf_func(i_s - 0.5, chord_positions[i_f])
            fil_wing.filaments_ym[i_s][i_f].start_coord = coord
            fil_wing.filaments_yp[i_s - 1][i_f].end_coord = coord
            if i_f == 1
                fil_wing.filaments_cw[i_s * 2 - 1][i_f].start_coord = coord
            elseif i_f < n_fil
                fil_wing.filaments_cw[i_s * 2 - 1][i_f].start_coord = coord
                fil_wing.filaments_cw[i_s * 2 - 1][i_f - 1].end_coord = coord
            else
                fil_wing.filaments_cw[i_s * 2 - 1][i_f - 1].end_coord = coord
            end
        end
    end
    for i in 1 : n_fil # tip y_plus
        coord = surf_func(n_strips + 1, chord_positions[i])
        fil_wing.filaments_yp[end][i].end_coord = coord
        if i == 1
            fil_wing.filaments_cw[end][i].start_coord = coord
        elseif i < n_fil
            fil_wing.filaments_cw[end][i].start_coord = coord
            fil_wing.filaments_cw[end][i - 1].end_coord = coord
        else
            fil_wing.filaments_cw[end][i - 1].end_coord = coord
        end
    end
    for i_f in 1 : n_fil # on the strip lines
        for i_s in 1 : n_strips
            coord = surf_func(i_s, chord_positions[i_f])
            fil_wing.filaments_ym[i_s][i_f].end_coord = coord
            fil_wing.filaments_yp[i_s][i_f].start_coord = coord
            if i_f == 1
                fil_wing.filaments_cw[i_s * 2][i_f].start_coord = coord
            elseif i_f < n_fil
                fil_wing.filaments_cw[i_s * 2][i_f].start_coord = coord
                fil_wing.filaments_cw[i_s * 2][i_f - 1].end_coord = coord
            else
                fil_wing.filaments_cw[i_s * 2][i_f - 1].end_coord = coord
            end
        end
    end
    return fil_wing
end

#= END VortexParticleWakeLUATATSolution functions =============================#

function transform_ThreeDVectors(
	func :: Function,
	coord :: ThreeDVector
	)
	return func(coord)
end

function transform_ThreeDVectors(
	func :: Function,
	fil :: ThreeDStraightVortexFilament
	)
	f = deepcopy(fil)
	f.start_coord = transform_ThreeDVectors(func, f.start_coord)
	f.end_coord = transform_ThreeDVectors(func, f.end_coord)
	return f
end

function transform_ThreeDVectors(
	func :: Function,
	coords :: Vector
	)
	return map(x->transform_ThreeDVectors(func, x), coords)
end

function transform_ThreeDVectors(
	func :: Function,
	fil_wing :: ThreeDSpanwiseFilamentWingRepresentation
	)
	w = deepcopy(fil_wing)
	for i = 1 : length(w.filaments_ym)
		w.filaments_yp[i] = transform_ThreeDVectors(func, w.filaments_yp[i])
		w.filaments_ym[i] = transform_ThreeDVectors(func, w.filaments_ym[i])
	end
    for i = 1 : length(w.filaments_cw)
        w.filaments_cw[i] = transform_ThreeDVectors(func, w.filaments_cw[i])
    end
	return w
end

function transform_ThreeDVectors(
	func :: Function,
	chord :: WingChordSection
	)
	c = deepcopy(chord)
	c.LE_location = transform_ThreeDVectors(func, c.LE_location)
	c.TE_location = transform_ThreeDVectors(func, c.TE_location)
	return c
end

function transform_ThreeDVectors(
	func :: Function,
	wing :: StripDefinedWing
	)
	w = deepcopy(wing)
	w.strips = map(x->transform_ThreeDVectors(func, x), w.strips)
	w.tip_yplus_LE_location =
		transform_ThreeDVectors(func, w.tip_yplus_LE_location)
	w.tip_yplus_TE_location =
		transform_ThreeDVectors(func, w.tip_yplus_TE_location)
	w.tip_yminus_LE_location =
		transform_ThreeDVectors(func, w.tip_yminus_LE_location)
	w.tip_yminus_TE_location =
		transform_ThreeDVectors(func, w.tip_yminus_TE_location)
	return w
end

#===============================================================================
    Integral remaps

    HJAB 2018
===============================================================================#

function telles_quadratic_remap(
    points :: Vector{Float64},
    weights :: Vector{Float64},
    singularity_location :: Float64)
    @assert(all(-1 .<= points .<= 1))
    @assert(abs(singularity_location) == 1)
    np =  (1 - points.*points) .* (singularity_location +
        sqrt(singularity_location*singularity_location - 1)) / 2 + points
    nw =  (-points .* (singularity_location +
        sqrt(singularity_location*singularity_location - 1)) + 1) .* weights
    return np, nw
end

function linear_remap(
    points :: Vector{T1},
    weights :: Vector{T1},
    original_range :: Vector{T2},
    new_range :: Vector{T3}
    ) where {T1 <: Real, T2 <: Number,  T3 <: Number}

    @assert(all(minimum(original_range) .<= points .<= maximum(original_range)))
    @assert(length(original_range) == 2)
    @assert(length(new_range) == 2)
    @assert(all(isfinite.(original_range)))
    @assert(all(isfinite.(new_range)))

    dx_ratio = (new_range[2] - new_range[1]) /
            (original_range[2] - original_range[1])
    offset = new_range[1] - original_range[1]
    nw = weights * dx_ratio
    np = (points - original_range[1]) * dx_ratio + offset + original_range[1]
    return np, nw
end

function exception_or_nonfinite_to_false(f::Function, arg::Float64)
    b = true
    try
        b = isfinite(f(arg))
    catch
        b = false
    end
    return b
end
#= END Integral remaps functions ==============================================#

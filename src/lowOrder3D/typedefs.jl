type KinemDef3D
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef

    function KinemDef3D(alpha :: MotionDef, h::MotionDef, u::MotionDef)
        new(alpha, h, u )
    end
end

immutable ThreeDFieldSimple
    f2d :: Vector{TwoDFlowField}
    function ThreeDFieldSimple()
        f2d = TwoDFlowField[]
        new(f2d)
    end
end

immutable ThreeDSurfSimple
    cref :: Float64
    AR :: Float64
    uref :: Float64
    pvt :: Float64
    lespcrit :: Vector{Float64}
    coord_file :: String
    ndiv :: Int8
    nspan :: Int8
    naterm :: Int8
    kindef :: KinemDef3D
    psi :: Vector{Float64}
    yle :: Vector{Float64}
    s2d :: Vector{TwoDSurf}
    a03d :: Vector{Float64}
    bc :: Vector{Float64}
    nshed :: Vector{Float64}
    bcoeff :: Vector{Float64}
    levstr :: Vector{Float64}
    fc :: Array{Float64}
    aterm3d :: Array{Float64}

    function ThreeDSurfSimple(AR, kindef, coord_file, pvt, lespcrit = [10.;]; nspan = 10, cref = 1., uref=1., ndiv=70, naterm=35)

        bref = AR*cref

        psi = zeros(nspan)
        yle = zeros(nspan)

        s2d = TwoDSurf[]

        for i = 1:nspan
            psi[i] = real(i)*(pi/2)/nspan
            yle[i] = -bref*cos(psi[i])/2.
        end

        #This code should be made more general to allow more motion types and combinations
        if typeof(kindef.h) == BendingDef
            for i = 1:nspan
                h_amp = evaluate(kindef.h.spl, yle[i])*kindef.h.scale
                h2d = CosDef(0., h_amp, kindef.h.k, kindef.h.phi)
                kinem2d = KinemDef(kindef.alpha, h2d, kindef.u)
                push!(s2d, TwoDSurf(coord_file, pvt,  kinem2d, lespcrit, c=cref, uref=uref, ndiv=ndiv, naterm=naterm))
            end
        else
            for i = 1:nspan
                kinem2d = KinemDef(kindef.alpha, kindef.h, kindef.u)
                lespc = lespcrit[1]
                push!(s2d, TwoDSurf(coord_file, pvt,  kinem2d, [lespc;], c=cref, uref=uref, ndiv=ndiv, naterm=naterm))
            end
        end

        a03d = zeros(nspan)
        aterm3d = zeros(naterm, nspan)

        bc = zeros(nspan)
        nshed = [0.;]
        bcoeff = zeros(nspan)
        levstr = zeros(nspan)
        fc = zeros(nspan,3)

        new(cref, AR, uref, pvt, lespcrit, coord_file,  ndiv, nspan, naterm, kindef,
        psi, yle, s2d, a03d, bc, nshed, bcoeff, levstr, fc, aterm3d)

    end
end

immutable KelvinConditionLLT
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
end

function (kelv::KelvinConditionLLT)(tev_iter::Array{Float64})
    val = zeros(kelv.surf.nspan)

    #Assume symmetry condition for now
    for i = 1:kelv.surf.nspan
        kelv.field.f2d[i].tev[end].s = tev_iter[i]

        #Update incduced velocities on airfoil
        update_indbound(kelv.surf.s2d[i], kelv.field.f2d[i])

        #Calculate downwash
        update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf.s2d[i])

        kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
        end

    calc_a0a13d(kelv.surf)

    for i = 1:kelv.surf.nspan
        val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i]
        + kelv.surf.a03d[i]) + 0.5*kelv.surf.aterm3d[1,i]

        nlev = length(kelv.field.f2d[i].lev)
        ntev = length(kelv.field.f2d[i].tev)

        for iv = 1:ntev
            val[i] = val[i] + kelv.field.f2d[i].tev[iv].s
        end
        for iv = 1:nlev
            val[i] = val[i] + kelv.field.f2d[i].lev[iv].s
        end
    end

    return val
end

immutable KelvinKuttaLLT
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
    nshed :: Int
end

function (kelv::KelvinKuttaLLT)(tev_iter::Array{Float64})
    val = zeros(kelv.surf.nspan + kelv.nshed)

    #Assume symmetry condition for now
    for i = 1:kelv.surf.nspan
        kelv.field.f2d[i].tev[end].s = tev_iter[i]
    end

    cntr = kelv.surf.nspan + 1
    for i = 1:kelv.surf.nspan
        if kelv.surf.s2d[i].levflag == 1
            kelv.field.f2d[i].lev[end].s = tev_iter[cntr]
            cntr += 1
        end
    end
    for i = 1:kelv.surf.nspan
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf.s2d[i], kelv.field.f2d[i])

        #Calculate downwash
        update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf.s2d[i])

        kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
    end

    calc_a0a13d(kelv.surf)



    for i = 1:kelv.surf.nspan
        val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i]
        + kelv.surf.a03d[i]) + 0.5*kelv.surf.aterm3d[1,i]

        nlev = length(kelv.field.f2d[i].lev)
        ntev = length(kelv.field.f2d[i].tev)

        for iv = 1:ntev
            val[i] = val[i] + kelv.field.f2d[i].tev[iv].s
        end
        for iv = 1:nlev
            val[i] = val[i] + kelv.field.f2d[i].lev[iv].s
        end
    end

    cntr = kelv.surf.nspan + 1
    for i = 1:kelv.surf.nspan
        if kelv.surf.s2d[i].levflag == 1
            if kelv.surf.s2d[i].a0[1] > 0
                lesp_cond = kelv.surf.s2d[i].lespcrit[1]
            else
                lesp_cond = -kelv.surf.s2d[i].lespcrit[1]
            end
            val[cntr] = kelv.surf.s2d[i].a0[1] + kelv.surf.a03d[1] - lesp_cond
            cntr += 1
        end
    end

    return val
end

#===============================================================================
    ThreeDVector

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDVector
    x :: Float64
    y :: Float64
    z :: Float64

    function ThreeDVector(x_ :: T1, y_ :: T2, z_ :: T3) where
        {T1 <: Real, T2 <: Real,  T3 <: Real}
        new(x_, y_, z_)
    end

    function ThreeDVector(array :: Vector{T}) where T <: Real
        @assert(size(array)[1] == 3)
        x = Float64(array[1])
        y = Float64(array[2])
        z = Float64(array[3])
        new(x, y, z)
    end
end

""" Convert type ThreeDVector to Array{Float64}"""
function convert(::Type{Array{T, 1}}, a::ThreeDVector) where T <: Real
    return [a.x, a.y, a.z]
end

""" Convert type Array{Real} to ThreeDVector"""
function convert(::Type{ThreeDVector}, a::Array{T, 1}) where T<: Real
    @assert(size(a)[1] == 3)
    b = ThreeDVector(a[1], a[2], a[3])
    return b
end

function Base.:+(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector([
        a.x + b.x,
        a.y + b.y,
        a.z + b.z ])
    return c
end

function Base.:-(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector([a.x - b.x, a.y - b.y, a.z - b.z])
    return c
end

function Base.:*(a::ThreeDVector, b::T) where T <: Real
    c = ThreeDVector([
        a.x * b,
        a.y * b,
        a.z * b ])
    return c
end

function Base.:*(a::T, b::ThreeDVector) where T <: Real
    return b * a
end

function Base.:/(a::ThreeDVector, b::T) where T <: Real
    c = ThreeDVector([
        a.x / b,
        a.y / b,
        a.z / b ])
    return c
end

"""Cross product of two ThreeDVectors"""
function cross(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector(
        [a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x ])
    return c
end

"""Dot product of two ThreeDVectors"""
function dot(a::ThreeDVector, b::ThreeDVector)
    return a.x * b.x + a.y * b.y + a.z * b.z
end

function Base.abs(a::ThreeDVector)
    return sqrt(a.x^2 + a.y^2 + a.z^2)
end

""" Set vector a to zero"""
function zero!(a::ThreeDVector)
    a.x = 0.0
    a.y = 0.0
    a.z = 0.0
    return Void
end

"""Return vector of length 1 with same direction"""
function unit(a::ThreeDVector)
    b = abs(a)
    return a / b
end

"""Set a vector normalised to length 1 with same direction"""
function unit!(a::ThreeDVector)
    b = abs(a)
    a.x /= b
    a.y /= b
    a.z /= b
    return Void
end

function Base.getindex(a::ThreeDVector, i::Int)  where Int <: Integer
    if i == 1 return a.x
    elseif i == 2 return a.y
    elseif i == 3 return a.z
    else throw(Core.BoundsError)
    end
end

function Base.size(a::ThreeDVector)
    return [3]
end

function Base.isfinite(a::ThreeDVector)
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z)
end
#= END ThreeDVector ----------------------------------------------------------=#

#===============================================================================
    ThreeDVortexParticle

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDVortexParticle
    coord :: ThreeDVector
    vorticity :: ThreeDVector
    size :: Float64

    velocity :: ThreeDVector
    vorticity_time_derivative :: ThreeDVector
end

"""
According to input vortex particle field particles, and kernal functions g and
f, compute the change for a timestep dt and return this as a new
particles_updated array of particles.
"""
function one_step(
    particles::Array{ThreeDVortexParticle},
    dt::Float64,
    g_function::Function,
    f_function::Function
    )

    delta_x = Array{ThreeDVector, 1}(size(particles))
    delta_vort = Array{ThreeDVector, 1}(size(particles))
    particles_updated = deepcopy(particles)

    function get_dx(particle_idx::Int64, particles::Array{ThreeDVortexParticle})
        location = particles[particle_idx].coord
        v = map(x->induced_velocity(x, location, g_function),
                        particles[vcat(1:particle_idx-1, particle_idx+1:end)])
        dx = sum(v) * dt
        return dx
    end
    function get_dvort(particle_idx::Int64, particles::Array{ThreeDVortexParticle})
        location = particles[particle_idx].coord
        vo = map(x->rate_of_change_of_vorticity(x, particles[particle_idx], g_function, f_function),
                    particles[vcat(1:particle_idx-1, particle_idx+1:end)])
        dvo = sum(vo) * dt
        return dvo
    end

    delta_x = map(i->get_dx(i, particles), 1:size(particles)[1])
    delta_vort = map(i->get_dvort(i, particles), 1:size(particles)[1])

    for i = 1 : size(particles)[1]
        particles_updated[i].coord += delta_x[i]
        particles_updated[i].vorticity += delta_vort[i]
    end
    return particles_updated
end

"""
Compute the velocity induced by ThreeDParticle particle at ThreeDCoordinate
coordinate, given a g function.
"""
function induced_velocity(
    particle::ThreeDVortexParticle,
    coordinate::ThreeDVector,
    g_function::Function
     )
     # Robertson 2010 Eq. 3
     rad = particle.coord - coordinate
     a = g_function(abs(rad)/particle.size) / (4. * pi)
     den = abs(rad)^3
     c = cross(rad, particle.vorticity)
     vel :: ThreeDVector = c * a / den
     return vel
end

"""
The vorticity of a particle over time changes as vortices stretch
 and what not. This RETURNS (doesn't change the value of) the rate of change
 of particle j with respect to time as domega_x / dt, domega_y / dt,
 domega_z / dt
"""
function rate_of_change_of_vorticity(
    particle_j::ThreeDVortexParticle,
    particle_k::ThreeDVortexParticle,
    g_function::Function,
    f_function::Function
    )
    # Robertson 2010 Eq. 4
    sigma_k = particle_k.size
    om_j = particle_j.vorticity
    om_k = particle_k.vorticity
    r = particle_j.coord - particle_k.coord
    g = g_function(abs(r) / particle_k.size)
    f = f_function(abs(r) / particle_k.size)
    om_cross = cross(particle_j.vorticity, particle_k.vorticity)
    r_cross_k = cross(r, particle_k.vorticity)
    om_dot_r = dot(particle_j.vorticity, r)

    t1 = 1. / (4 * pi * sigma_k ^ 3)
    t21n = - g * om_cross
    t21d = (abs(r) / sigma_k) ^ 3
    t21 = t21n / t21d
    t221 = 1 / (abs(r)^2)
    t222 = 3 * g / (abs(r) / sigma_k) ^ 3 - f
    t223 = om_dot_r * r_cross_k
    t22 = t221 * t222 * t223
    t2 = t21 + t22
    t = t1 * t2
    return t
end

"""Take a set of vortex particles and the value of g for each of them,
 and computes the velocity of each. Returns a vector of velocities. """
function particle_velocities(
    particles::Array{ThreeDVortexParticle},
    particle_g_function::Function)

    @assert(particles.size() == particles_g_values.size())  # Same size arrays

    vel = zeros(particles.size())
    i :: Int64
    v0 :: ThreeDVector
    for i = 1 : particles.size()
        x0 = particles[i].coord
        zero!(v0)
        vel[i] = mapreduce(x->induced_velocity(x, x0, particle_g_function),
            +, v0, particles)
    end
    return vel
end

"""Take a set of vortex particles and the value of g and f for each of them,
 and computes the rate of change of vorticity for each of them. """
function particle_rate_of_change_of_vorticity(
    particles::Array{ThreeDVortexParticle},
    particle_g_function::Function,
    particle_f_function::Function)

    @assert(particles.size() == particles_g_values.size())  # Same size arrays

    dvort = zeros(particles.size())
    i :: Int64
    v0 :: ThreeDVector
    for i = 1 : particles.size()
        x0 = particles[i].coord
        zero!(v0)
        dvort[i] = mapreduce(
            x->rate_of_change_of_vorticity(dvort[i], x[1],
                                        particle_g_values, particle_f_values),
            +, v0, particles)
    end
    return vel
end
#= END ThreeDVortexParticle --------------------------------------------------=#

#===============================================================================
    ThreeDStraightVortexFilament

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDStraightVortexFilament
    start_coord :: ThreeDVector
    end_coord :: ThreeDVector
    vorticity :: Float64

    function ThreeDStraightVortexFilament(
        start_coord :: ThreeDVector,
        end_coord :: ThreeDVector,
        vorticity :: T
        ) where T <: Real
        new(start_coord, end_coord, Float64(vorticity))
    end
end

function ThreeDStraightVortexFilament(
    start_coord :: ThreeDVector,
    end_coord :: ThreeDVector
    )
    return ThreeDStraightVortexFilament(start_coord, end_coord, 0.0)
end

ThreeDStraightVortexFilament() =
    ThreeDStraightVortexFilament(
        ThreeDVector([0., 0., 0.]),
        ThreeDVector([1., 1., 1.]),
        0.0)
#= END ThreeDStraightVortexFilament ------------------------------------------=#

#===============================================================================
    ThreeDVortexParticleSet

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDVortexParticleSet
    particles :: Vector{ThreeDVortexParticle}

    _reduction_factor_fn :: Function
    _vorticity_fraction_fn :: Function
end
#= END ThreeDVortexParticleSet -----------------------------------------------=#


#===============================================================================
    WingChordSection

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type WingChordSection
    # Define the locations of the strips used in the analysis:
    LE_location :: ThreeDVector
    TE_location :: ThreeDVector
    camber_line :: Spline1D # Y locations for x in [-1, 1]

    function WingChordSection(
        LE_location :: ThreeDVector,
        TE_location :: ThreeDVector,
        camber_line :: Spline1D
        )
        new(LE_location, TE_location, camber_line)
    end
end

function WingChordSection(
    LE_location :: ThreeDVector,
    TE_location :: ThreeDVector
    )
    xs = [-1., 1.]
    ys = [0., 0.]
    camber_line = Spline1D(xs, ys, k=1)
    return WingChordSection(LE_location, TE_location, camber_line)
end

function WingChordSection()
    return WingChordSection(
        ThreeDVector(-1., 0., 0.),
        ThreeDVector(1., 0., 0.)
    )
end

function angle_of_attack(chord :: WingChordSection)
    dx = chord.TE_location - chord.LE_location
    unit!(dx)
    math.asin(dx.z)
end

"""Compute the location for the camberline at a point x in [-1, 1] for [le, te]
"""
function location(
    chord_section :: WingChordSection,
    normal_dir :: ThreeDVector,
    x :: Float64)

    dx = chord_section.TE_location - chord_section.LE_location
    len = abs(dx)
    cam = len * chord_section.camber_line(x) / 2
    loc = chord_section.LE_location + dx * (x + 1) / 2
    loc += unit(normal_dir) * cam
    return loc
end

function chord(
    chord_section :: WingChordSection
    )

    return chord_section.TE_location - chord_section.LE_location
end

function lump_vorticities(
    chord_section :: WingChordSection,
    vorticity_fn :: Function,
    locations :: Vector{Float64},
    extenal_effects :: Vector{Float64}
    )
    @assert(all(map(x->abs(x) <= 1.0, locations)))
    @assert(allunique(locations))
    @assert(size(locations) == size(extenal_effects))
    # We want things in order, but need to return our results in whatever order
    # the user needs:
    reordering = sortperm(locations)
    sorted = locations[reordering]
    separators = vcat(-1.0, (sorted[1 : end - 1] + sorted[2 : end]) / 2., 1.0)
    values = zeros(size(locations)[1])
    # We need to correct our integral for chord length and for the camber
    len = abs(chord(chord_section))
    function fn(x :: Float64)
        mult = len * sqrt(1 + (derivative(chord_section.camber_line, x))^2)
        return mult * vorticity_fn(x)
    end
    # Use the Simpson's rule over each separated bit to integrate a vorticity.
    for i = 1 : size(locations)[1]
        x = [separators[i], 0.0, separators[i + 1]]
        x[2] = (x[1] + x[3])/ 2.
        w = (x[3] - x[1]) * [1.0, 4.0, 1.0] / 3.0
        values[reordering[i]] = extenal_effects[i] *
            mapreduce(x->x[1] * fn(x[2]), +, 0.0, zip(w, x))
    end
    return values
end

function lump_vorticities(
    chord_section :: WingChordSection,
    vorticity_fn :: Function,
    locations :: Vector{Float64}
    )
    ext = ones(size(locations)[1])
    return lump_vorticities(chord_section, vorticity_fn, locations, ext)
end
#= END WingChordSection ------------------------------------------------------=#

#===============================================================================
    StripDefinedWing

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type StripDefinedWing
    # Define the locations of the strips used in the analysis:
    strips :: Vector{WingChordSection}

    # Define the locations of the wing tips:
    tip_yplus_LE_location :: ThreeDVector
    tip_yplus_TE_location :: ThreeDVector
    tip_yminus_LE_location :: ThreeDVector
    tip_yminus_TE_location :: ThreeDVector
end

""" Obtain 3 splines representing the x, y and z positions respectively
of the wing trailing edge with respect to their index. For n segements, 0
represents the yminus tip and n + 1 the yplus tip. """
function te_spline(wing :: StripDefinedWing)
    n = size(wing.strips)[1]
    x = zeros(n + 2)
    y = zeros(n + 2)
    z = zeros(n + 2)
    x[1] = wing.tip_yminus_TE_location.x
    y[1] = wing.tip_yminus_TE_location.y
    z[1] = wing.tip_yminus_TE_location.z
    for i = 1:n
        if isfinite(wing.strips[i].TE_location) != true
            error("Nonfinite wing strip TE definition.")
        end
        x[i+1] = wing.strips[i].TE_location.x
        y[i+1] = wing.strips[i].TE_location.y
        z[i+1] = wing.strips[i].TE_location.z
    end
    x[n+2] = wing.tip_yplus_TE_location.x
    y[n+2] = wing.tip_yplus_TE_location.y
    z[n+2] = wing.tip_yplus_TE_location.z
    iota_array = Vector{Float64}(0:n+1)
    spl_x = Spline1D(iota_array, x)
    spl_y = Spline1D(iota_array, y)
    spl_z = Spline1D(iota_array, z)
    return spl_x, spl_y, spl_z
end

""" Obtain 3 splines representing the x, y and z positions respectively
of the wing leading edge with respect to their index. For n segements, 0
represents the yminus tip and n + 1 the yplus tip. """
function le_spline(wing :: StripDefinedWing)
    n = size(wing.strips)[1]
    x = zeros(n + 2)
    y = zeros(n + 2)
    z = zeros(n + 2)
    x[1] = wing.tip_yminus_LE_location.x
    y[1] = wing.tip_yminus_LE_location.y
    z[1] = wing.tip_yminus_LE_location.z
    for i = 1:n
        if isfinite(wing.strips[i].LE_location) != true
                error("Nonfinite wing strip LE definition.")
            end
        x[i+1] = wing.strips[i].LE_location.x
        y[i+1] = wing.strips[i].LE_location.y
        z[i+1] = wing.strips[i].LE_location.z
    end
    x[n+2] = wing.tip_yplus_LE_location.x
    y[n+2] = wing.tip_yplus_LE_location.y
    z[n+2] = wing.tip_yplus_LE_location.z
    iota_array = Vector{Float64}(0:n+1)
    spl_x = Spline1D(iota_array, x)
    spl_y = Spline1D(iota_array, y)
    spl_z = Spline1D(iota_array, z)
    return spl_x, spl_y, spl_z
end

""" Obtain the direction of the normal vector excluding any camber.
Returns 3 splines representing the x, y and z direction in a unit vector
representing the normal a the midchord for arg is 0 to n+1 for n strips
on the wing. """
function nocamber_normal_splines(wing :: StripDefinedWing)
    lex, ley, lez = le_spline(wing)
    tex, tey, tez = te_spline(wing)
    n = size(wing.strips)[1]
    x = zeros(n + 2)
    y = zeros(n + 2)
    z = zeros(n + 2)
    for i = 0 : n + 1
        le = ThreeDVector(lex(i), ley(i), lez(i))
        te = ThreeDVector(tex(i), tey(i), tez(i))
        dle = ThreeDVector(map(x->derivative(x, Float64(i)), [lex, ley, lez]))
        dte = ThreeDVector(map(x->derivative(x, Float64(i)), [tex, tey, tez]))
        cdir = te - le # chord direction
        if abs(cdir) == 0.0
            normal = ThreeDVector(0.0, 0.0, 0.0)
        else
            ddir = (dle + dte) / 2. # midchord spanwise direction
            normal = unit(cross(cdir, ddir))
        end
        x[i + 1] = normal.x
        y[i + 1] = normal.y
        z[i + 1] = normal.z
    end
    iota_array = Vector{Float64}(0 : n+1)
    spl_x = Spline1D(iota_array, x)
    spl_y = Spline1D(iota_array, y)
    spl_z = Spline1D(iota_array, z)
    return spl_x, spl_y, spl_z
end

"""
Returns a function that generates points on the wing surface.

For a StripDefinedWing it is useful to be able to obtain a continious surface.
This function returns a function f(s, x) that returns a point p (ThreeDVector)
on the wing surface. s is the strip position (0 -> y_minus tip, n + 1 to y_plus
tip) and x defines the chordwise position.
"""
function get_surface_fn(
    wing :: StripDefinedWing
    )
    lex, ley, lez = le_spline(wing)
    tex, tey, tez = te_spline(wing)
    nex, ney, nez = nocamber_normal_splines(wing)

    function s(
        strip_pos :: T1,
        x :: T2
        ) where {T1 <: Real, T2 <: Real}
        @assert(abs(x) <= 1.)
        @assert(strip_pos >= 0)
        @assert(strip_pos <= size(wing.strips)[1] + 1)
        const s = strip_pos
        le = ThreeDVector(lex(s), ley(s), lez(s))
        te = ThreeDVector(tex(s), tey(s), tez(s))
        n = ThreeDVector(nex(s), ney(s), nez(s))
        p = le + 0.5 * (x + 1.) * (te - le)
        i_cf = Int64(floor(strip_pos))
        i_cc = Int64(ceil(strip_pos))
        if 0 < i_cf <= length(wing.strips)
            c_cf = wing.strips[i_cf].camber_line(x)
        else
            c_cf = 0.0
        end
        if 0 < i_cc <= length(wing.strips)
            c_cc = wing.strips[i_cc].camber_line(x)
        else
            c_cc = 0.0
        end
        p += n * (c_cc * (x % 1.0) + c_cf * (1. - x % 1.0)) * abs(te - le)
        if isfinite(p) != true
            error("Evaluated surface location as non-finite")
        end
        return p
    end
    return s
end

"""
Returns a function that returns the direction of the wing chord
"""
function get_chord_dir_fn(
    wing :: StripDefinedWing
    )
    lex, ley, lez = le_spline(wing)
    tex, tey, tez = te_spline(wing)
    function s(
        strip_pos :: T1,
        x :: T2
        ) where {T1 <: Real, T2 <: Real}
        @assert(abs(x) <= 1.)
        @assert(strip_pos >= 0)
        @assert(strip_pos <= size(wing.strips)[1] + 1)
        const s = strip_pos
        le = ThreeDVector(lex(s), ley(s), lez(s))
        te = ThreeDVector(tex(s), tey(s), tez(s))
        cdir = unit(te - le)
        if isfinite(cdir) != true
            error("Tried to evaluate chord direction at zero-chord location.")
        end
        return p
    end
    return s
end

"""
Returns a function that returns the deta_dx dot unit(chord)
"""
function get_surface_detadx_dot_c_fn(
    wing :: StripDefinedWing
    )
    lex, ley, lez = le_spline(wing)
    tex, tey, tez = te_spline(wing)
    nex, ney, nez = nocamber_normal_splines(wing)

    function s(
        strip_pos :: T1,
        x :: T2
        ) where {T1 <: Real, T2 <: Real}
        @assert(abs(x) <= 1.)
        @assert(strip_pos >= 0)
        @assert(strip_pos <= size(wing.strips)[1] + 1)
        const s = strip_pos
        le = ThreeDVector(lex(s), ley(s), lez(s))
        te = ThreeDVector(tex(s), tey(s), tez(s))
        n = ThreeDVector(nex(s), ney(s), nez(s))
        cdir = unit(te - le)
        i_cf = Int64(floor(strip_pos))
        i_cc = Int64(ceil(strip_pos))
        if 0 < i_cf <= length(wing.strips)
            c_cf = derivative(wing.strips[i_cf].camber_line(x))
        else
            c_cf = 0.0
        end
        if 0 < i_cc <= length(wing.strips)
            c_cc = derivative(wing.strips[i_cc].camber_line(x))
        else
            c_cc = 0.0
        end
        cderiv = c_cc * (x % 1.0) + c_cf * (1. - x % 1.0)
        deta_dx = cderiv * abs(te -le) / 2.0
        p = deta_dx * cdir
        if isfinite(p) != true
            error("Could not evaluate deta/dx * unit(chord)")
        end
        return p
    end
    return s
end

#= END StripDefinedWing ------------------------------------------------------=#

#===============================================================================
    ThreeDSpanwiseFilamentWingRepresentation

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDSpanwiseFilamentWingRepresentation
    filaments_ym :: Vector{Vector{ThreeDStraightVortexFilament}}
    filaments_yp :: Vector{Vector{ThreeDStraightVortexFilament}}
    n_filaments_per_chord :: Vector{Int64}
    n_chords :: Int64

    function ThreeDSpanwiseFilamentWingRepresentation(
        n_chords :: Int64, n_fils_per_chord :: Int64
        )

        n_filaments_per_chord = ones(n_chords) * n_fils_per_chord
        c_fils = Vector{ThreeDStraightVortexFilament}(n_fils_per_chord)
        for i = 1 : length(c_fils)
            c_fils[i] = ThreeDStraightVortexFilament()
        end
        filaments_ym =
            [deepcopy(c_fils) for _ in 1:n_chords]
        filaments_yp = deepcopy(filaments_ym)
        new(filaments_ym, filaments_yp, n_filaments_per_chord, n_chords)
    end
end

function convert(
    ::Type{Vector{ThreeDStraightVortexFilament}},
    a::ThreeDSpanwiseFilamentWingRepresentation)

    n_fil = 2 * sum(a.n_filaments_per_chord)
    vect = Vector{ThreeDStraightVortexFilament}([])
    for i = 1 : n_chords
        vect = vcat(vect, a.filaments_ym[i], a.filaments_yp[i])
    end
    return vect
end

function zero_vorticities!(wing :: ThreeDSpanwiseFilamentWingRepresentation)
    for c in vcat(filaments_ym, filaments_yp)
        for f in c
            c.vorticity = 0.0
        end
    end
end

function add_vorticity!(
    wing :: StripDefinedWing,
    filament_positions :: Vector{Float64},
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    strip_idx :: Int64,
    func :: Function)
    @assert(strip_idx > 0)
    @assert(strip_idx <= length(wing.strips))
    @assert(length(fil_wing.filaments_yp[strip_idx]) ==
        length(fil_wing.filaments_ym[strip_idx]) )

    nf = length(fil_wing.filaments_ym[strip_idx])
    fils_yp = fil_wing.filaments_yp[strip_idx]
    fils_ym = fil_wing.filaments_ym[strip_idx]

    ext_ym, ext_yp = slant_correction_factors(
        wing, fil_wing, filament_positions, strip_idx)
    vort = lump_vorticities(wing.strips[strip_idx], func,
        filament_positions, ext_ym)
    for i = 1 : nf
        fil_wing.filaments_ym[strip_idx][i].vorticity = vort[i]
    end
    vort = lump_vorticities(wing.strips[strip_idx], func,
        filament_positions, ext_yp)
    for i = 1 : nf
        fil_wing.filaments_yp[strip_idx][i].vorticity = vort[i]
    end
    return
end
#= END ThreeDSpanwiseFilamentWingRepresentation ------------------------------=#

#===============================================================================
    VortexParticleWakeLAUTATSolution

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type VortexParticleWakeLAUTATSolution
    wing :: StripDefinedWing
    wake :: ThreeDVortexParticleSet
    filament_wing :: ThreeDSpanwiseFilamentWingRepresentation

    free_stream_velocity :: ThreeDVector
    # external_purturbation: accepts ThreeDVector coord & returns ThreeDVector.
    external_purturbation :: Function

    time :: Float64

    n_fourier_terms :: Int64
    fourier_terms :: Vector{Vector{Float64}}
    old_fourier_terms :: Vector{Vector{Float64}}

    k_sloc :: Vector{Float64}
    k_sind :: Vector{Int64}
    old_fil_wing_bound_vorticity_vector :: Vector{Float64}

    function VortexParticleWakeLAUTATSolution()
        wing = StripDefinedWing()
        wake = ThreeDVortexParticleSet()
        filament_wing = ThreeDSpanwiseFilamentWingRepresentation()
        free_stream_velocity = ThreeDVector(1, 0, 0)
        new( wing, wake, filament_wing, free_stream_velocity)
    end
end
#= END VortexParticleWakeLAUTATSolution --------------------------------------=#

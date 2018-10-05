mutable struct KinemDef3D
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef

    function KinemDef3D(alpha :: MotionDef, h::MotionDef, u::MotionDef)
        new(alpha, h, u )
    end
end

struct ThreeDFieldSimple
    f2d :: Vector{TwoDFlowField}
    function ThreeDFieldSimple()
        f2d = TwoDFlowField[]
        new(f2d)
    end
end

struct ThreeDSurfSimple
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

    function StripDefinedWing(
        strips :: Vector{WingChordSection},
        tip_yplus_LE_location :: ThreeDVector,
        tip_yplus_TE_location :: ThreeDVector,
        tip_yminus_LE_location :: ThreeDVector,
        tip_yminus_TE_location :: ThreeDVector,
        )

        new(strips, tip_yplus_LE_location, tip_yplus_TE_location,
            tip_yminus_LE_location, tip_yminus_TE_location)
    end
end

function StripDefinedWing()
    strips = [WingChordSection()]
    return StripDefinedWing(
        strips, ThreeDVector(-1, -1, 0), ThreeDVector(1, -1, 0),
        ThreeDVector(-1, 1, 0), ThreeDVector(1, 1, 0)
    )
end

""" Obtain 3 splines representing the x, y and z positions respectively
of the wing trailing edge with respect to their index. For n segements, 0.5
represents the yminus tip and n + .5 the yplus tip. """
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
    iota_array = vcat(0.5, 1:n, n+0.5)
    s_o = min(3, length(x) - 1)
    spl_x = Spline1D(iota_array, x, k = s_o)
    spl_y = Spline1D(iota_array, y, k = s_o)
    spl_z = Spline1D(iota_array, z, k = s_o)
    return spl_x, spl_y, spl_z
end

""" Obtain 3 splines representing the x, y and z positions respectively
of the wing leading edge with respect to their index. For n segements, 0.5
represents the yminus tip and n + 0.5 the yplus tip. """
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
    iota_array = vcat(0.5, 1:n, n+0.5)
    s_o = min(3, length(x) - 1)
    spl_x = Spline1D(iota_array, x, k = s_o)
    spl_y = Spline1D(iota_array, y, k = s_o)
    spl_z = Spline1D(iota_array, z, k = s_o)
    return spl_x, spl_y, spl_z
end

""" Obtain the direction of the normal vector excluding any camber.
Returns 3 splines representing the x, y and z direction in a unit vector
representing the normal a the midchord for arg is .5 to n+.5 for n strips
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
    iota_array = vcat(0.5, 1:n, n+0.5)
    s_o = min(3, length(x) - 1)
    spl_x = Spline1D(iota_array, x, k = s_o)
    spl_y = Spline1D(iota_array, y, k = s_o)
    spl_z = Spline1D(iota_array, z, k = s_o)
    return spl_x, spl_y, spl_z
end

"""
Returns a function that generates points on the wing surface.

For a StripDefinedWing it is useful to be able to obtain a continious surface.
This function returns a function f(s, x) that returns a point p (ThreeDVector)
on the wing surface. s is the strip position (0.5 -> y_minus tip, n + .5 to
y_plus tip) and x defines the chordwise position.
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
        @assert(strip_pos >= 0.5)
        @assert(strip_pos <= size(wing.strips)[1] + .5)
        const s = strip_pos
        le = ThreeDVector(lex(s), ley(s), lez(s))
        te = ThreeDVector(tex(s), tey(s), tez(s))
        cdir = unit(te - le)
        if isfinite(cdir) != true
            error("Tried to evaluate chord direction at zero-chord location.")
        end
        return cdir
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
            c_cf = derivative(wing.strips[i_cf].camber_line, x)
        else
            c_cf = 0.0
        end
        if 0 < i_cc <= length(wing.strips)
            c_cc = derivative(wing.strips[i_cc].camber_line, x)
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
    filaments_cw :: Vector{Vector{ThreeDStraightVortexFilament}}

    function ThreeDSpanwiseFilamentWingRepresentation(
        n_chords :: Int64, n_fils_per_chord :: Int64
        )

        n_filaments_per_chord = ones(n_chords) * n_fils_per_chord
        c_fils = Vector{ThreeDStraightVortexFilament}(n_fils_per_chord)
        for i = 1 : length(c_fils)
            c_fils[i] = ThreeDStraightVortexFilament()
        end
        filaments_ym = [deepcopy(c_fils) for _ in 1:n_chords]
        filaments_yp = deepcopy(filaments_ym)
        cw_fils = c_fils[1:end - 1]
        filaments_cw = [deepcopy(cw_fils) for _ in 1:n_chords * 2 + 1]
        new(filaments_ym, filaments_yp, filaments_cw)
    end
end

function convert(
    ::Type{Vector{ThreeDStraightVortexFilament}},
    a::ThreeDSpanwiseFilamentWingRepresentation)

    vect = Vector{ThreeDStraightVortexFilament}([])
    for i = 1 : length(a.filaments_ym)
        vect = vcat(vect, a.filaments_ym[i], a.filaments_yp[i])
    end
    for v in a.filaments_cw
        vect = vcat(vect, v)
    end
    return vect
end

function zero_vorticities!(wing :: ThreeDSpanwiseFilamentWingRepresentation)
    for c in vcat(wing.filaments_ym, wing.filaments_yp, wing.filaments_cw)
        for f in c
            f.vorticity = 0.0
        end
    end
end

function add_vorticity!(
    wing :: StripDefinedWing,
    filament_positions :: Vector{Float64},  # in  [-1,1]
    fil_wing :: ThreeDSpanwiseFilamentWingRepresentation,
    strip_idx :: Int64,
    func :: Function)   # For x in [-1, 1]
    @assert(strip_idx > 0)
    @assert(strip_idx <= length(fil_wing.filaments_yp))
    @assert(length(fil_wing.filaments_yp[strip_idx]) ==
        length(fil_wing.filaments_ym[strip_idx]) )
    @assert(all(-1 .< filament_positions .< 1))
    @assert(isfinite(func(0.0)))

    nf = length(fil_wing.filaments_ym[strip_idx])
    fils_yp = fil_wing.filaments_yp[strip_idx]
    fils_ym = fil_wing.filaments_ym[strip_idx]

    ext_ym, ext_yp = slant_correction_factors(
        wing, fil_wing, filament_positions, strip_idx)
    vort = lump_vorticities(wing.strips[strip_idx], func,
        filament_positions, ext_ym)
    for i = 1 : nf
        fil_wing.filaments_ym[strip_idx][i].vorticity += vort[i]
    end
    for i = 1 : nf - 1
        fil_wing.filaments_cw[strip_idx * 2 - 1][i].vorticity += sum(vort[1:i])
        fil_wing.filaments_cw[strip_idx * 2][i].vorticity -= sum(vort[1:i])
    end

    vort = lump_vorticities(wing.strips[strip_idx], func,
        filament_positions, ext_yp)
    for i = 1 : nf
        fil_wing.filaments_yp[strip_idx][i].vorticity += vort[i]
    end
    for i = 1 : nf - 1
        fil_wing.filaments_cw[strip_idx * 2 + 1][i].vorticity -= sum(vort[1:i])
        fil_wing.filaments_cw[strip_idx * 2][i].vorticity += sum(vort[1:i])
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

    # The shedding locations as a function of span position
    k_sloc :: Vector{Float64}
    # The index associated with each shedding location
    k_sind :: Vector{Int64}
    # The old bound vorticity vector
    old_fil_wing_bound_vorticity_vector :: Vector{Float64}

    function VortexParticleWakeLAUTATSolution()
        wing = StripDefinedWing()
        wake = ThreeDVortexParticleSet()
        filament_wing = build_vortex_filament_wing_geometry(wing, [-1., 0., 1])
        free_stream_velocity = ThreeDVector(1, 0, 0)
        new( wing, wake, filament_wing, free_stream_velocity)
    end
end
#= END VortexParticleWakeLAUTATSolution --------------------------------------=#

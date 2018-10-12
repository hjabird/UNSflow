#===============================================================================
    Lautat3D.jl

    A three dimensional version of the Ramesh Large amplitude unsteady
    thin-aerofoil theory.

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

mutable struct Lautat3D
    # Physical definitions:
    wing :: EquationSurf    # x chordwise (-1 = LE), y spanwise
    kinematics :: CoordinateTransform3D

    timestep :: Float64
    strip_centres :: Vector{Float64}
    strip_divisions :: Int64
    strip_discretisation :: Vector{Float64}
    num_fourier_coeffs :: Int64

    shed_particle_radius :: Float64
    shed_particle_kernal :: Vortex3DRegularisationFunctions

    # Problem solving state:
    wake :: Vorticity3DSimpleCollector
    last_shed_tevrs :: Vorticity3DSimpleCollector
    wing_geometry :: Vector{Vector{BilinearQuad}}
    wing_aero_discretisation :: Vorticity3DSimpleCollector

    # We divide the space to either side of each strip at least into 2.
    # Hence we need to store the location of all these divides + the wing tips.
    extended_strip_centers :: Vector{Float64}

    # We'll be the geometry data for each strip again and again, so we'll cache
    # it. These points are in [-1, 1], including the endpoints. It is easy to
    # compute the centres of the vortex rings from them.
    strip_eval_points :: Vector{Vector{Vector3D}}
    strip_eval_normals :: Vector{Vector{Vector3D}}
    strip_eval_points_velocities :: Vector{Vector{Vector3D}}

    strip_fourier_coeffs :: Vector{Vector{Float64}}
    old_strip_fourier_coeffs :: Vector{Vector{Float64}}

    function Lautat3D()
        m_wing = EquationSurf(x->Vector3D(x[1], x[2], 0))
        m_kinematics = CoordinateTransform3D()
        timestep = NaN
        m_strip_centres = Float64[]
        strip_divisions = 0
        strip_discretisation = []
        num_fourier_coeffs = 0
        shed_particle_radius = -1.
        shed_particle_kernal = threed_planetary_kernels()
        wake = Vorticity3DSimpleCollector()
        last_shed_tevrs = Vorticity3DSimpleCollector()
        wing_geometry = Vector{Vector{BilinearQuad}}()
        wing_aero_discretisation = Vorticity3DSimpleCollector()
        extended_strip_centers = Float64[]
        strip_eval_points = Vector{Vector{Vector3D}}()
        strip_eval_points_velocities = Vector{Vector{Vector3D}}()
        strip_eval_normals = Vector{Vector{Vector3D}}()
        strip_fourier_coeffs = Vector{Vector{Float64}}()
        old_strip_fourier_coeffs = Vector{Vector{Float64}}()
        new(m_wing, m_kinematics, timestep,
            m_strip_centres, strip_divisions, strip_discretisation,
            num_fourier_coeffs, shed_particle_radius,
            shed_particle_kernal, wake, last_shed_tevrs, wing_geometry,
            wing_aero_discretisation, extended_strip_centers,
            strip_eval_points, strip_eval_normals,
            strip_fourier_coeffs, old_strip_fourier_coeffs)
    end
end

#= DISCRETISATION ------------------------------------------------------------=#

# Compute a.wing_geometry
function discretise_wing_to_geometry!(a::Lautat3D)
    @assert(length(a.strip_discretisation) > 10, string(
        "length(Lautat3D.strip_discretisation) probably needs to be more than",
        " 10.", " Current value is ", length(a.strip_discretisation), "."))
    @assert(extrema(a.strip_discretisation) == (-1, 1), string("All of the ",
        "Lautat3D.strip_discretisation values should be -1 <= val <= 1. ",
        "There should be points at -1 and 1. ",
        "strip_discretisation was ", a.strip_discretisation))
    @assert(issorted(a.strip_discretisation), "Expected "*
        "Lautat3D.strip_discretisation to be sorted in ascending order.")

    generate_extended_strip_centres!(a)
    a.wing_geometry = Vector{Vector{BilinearQuad}}(undef,
        2 * a.strip_divisions * length(a.strip_centres))
    geom = discretise(a.wing, BilinearQuad, a.strip_discretisation,
        a.extended_strip_centers)
    for i = 1 : length(a.extended_strip_centers) - 1
        a.wing_geometry[i] = geom[:, i]
    end
    return
end

# Turn a.wing_geometry into an aerodynamic object.
function generate_wing_aero_lattice_from_geometry!(a::Lautat3D)
    ns = length(a.extended_strip_centers)
    @assert(ns - 1 == length(a.wing_geometry), string("Expected Lautat3D",
        " internals extended_strip_centers to be one longer than ",
        "wing geometry. Perhaps discretise_wing_to_geometry! hasn't",
        " been called? Lengths were ", ns, " and ", length(a.wing_geometry),
        "."))
    a.wing_aero_discretisation = Vorticity3DSimpleCollector()
    for i = 1 : length(a.wing_geometry)
        new_collector = Vorticity3DSimpleCollector()
        for j = 1 : length(a.wing_geometry[i])
            push!(new_collector, VortexRing(a.wing_geometry[i][j]))
        end
        push!(a.wing_aero_discretisation, new_collector)
    end
    return
end

function generate_extended_strip_centres!(a::Lautat3D)
    @assert(length(a.strip_centres) >= 1, string("length(a.strip_centres) = ",
        length(a.strip_centres), " but should be more than or equal to 1."))
    @assert(all(-1 .< a.strip_centres .< 1), string("All of the ",
        "Lautat3D.strip_centres values should be -1 < val < 1. ",
        "strip_centres was ", a.strip_centres))
    @assert(issorted(a.strip_centres), "Expected Lautat3D.strip_centres"*
        " to be sorted in ascending order.")
    @assert(a.strip_divisions > 0, string("Strip divisions must be more than ",
        "zero. Currently ", a.strip_divisions, "."))

    n_strips = 2 * a.strip_divisions * length(a.strip_centres)
    a.extended_strip_centers = Vector{Float64}(undef, n_strips + 1)
    for i = 0 : n_strips
        # We want a value that will interpolate a.strip_centres. So strip_pos
        # == 1 would be on the a.strip_centres[1], 1.5 half way between
        # stip_centres[1] and strip_centres[2].
        strip_pos = i / (2 * a.strip_divisions) + 0.5
        strip_upper = Int64(ceil(strip_pos))
        strip_lower = Int64(floor(strip_pos))
        strip_y_upper = strip_upper > length(a.strip_centres) + 0.5 ? 1. :
            a.strip_centres[strip_upper]
        strip_y_lower = strip_lower == 0 ? -1. : a.strip_centres[strip_lower]
        # Width of the bit we're interpolating.
        h = (strip_y_upper == 1) || (strip_y_lower == -1) ? 0.5 : 1.
        strip_pos = (strip_pos % 1) / h
        strip_pos = strip_y_lower == -1 ? strip_pos - 1 : strip_pos
        a.extended_strip_centers[i+1] = strip_pos * strip_y_upper +
            (1-strip_pos) * strip_y_lower
    end
    return
end

function generate_eval_points_data!(a::Lautat3D)
    evl_x = (a.strip_discretisation[1:end-1] .+a.strip_discretisation[2:end]) .*
        0.5
    tmp_wing = EquationSurf(x->evaluate(a.kinematics, evaluate(a.wing, x)))
    nsc = length(a.strip_centres)
    a.strip_eval_points = Vector{Vector{Vector3D}}(undef, nsc)
    a.strip_eval_normals = Vector{Vector{Vector3D}}(undef, nsc)
    a.strip_eval_points_velocities = Vector{Vector{Vector3D}}(undef, nsc)
    for i = 1 : nsc
        a.strip_eval_points[i] =
            map(x->evaluate(tmp_wing, [x, a.strip_centres[i]]), evl_x)
        a.strip_eval_normals[i] = map(x->unit(normal(tmp_wing,
            [x, a.strip_centres[i]])), evl_x)
        a.strip_eval_points_velocities[i] = map(x->derivative(a.kinematics,
            evaluate(a.wing, [x, a.strip_centres[i]])), evl_x)
    end
    return
end

function shed_new_tevrs(a::Lautat3D)
    @assert(length(a.old_strip_fourier_coeffs) ==
        length(a.strip_fourier_coeffs), "Non-matching fourier vector lengths")
    # New trailing edge vortex rings
    new_tevrs = Vector{VortexRing}(undef, length(a.strip_discretisation))
    old_tevrs = deepcopy(a.last_shed_tevrs)
    old_tevr_vorticities = vorticity_vector(old_tevs)
    if length(a.last_shed_tevrs) != length(new_tevrs)
        # We need to base our new geometry on induced_velocity.
        for i = 1 : length(new_tevrs)
            c1 = coords(a.wing_aero_discretisation[i][end].geometry)[2]
            c4 = coords(a.wing_aero_discretisation[i][end].geometry)[3]
            c2 = (induced_velocity(a.wing_aero_discretisation, c1) +
                induced_velocity(a.wake_discretisation, c1)) * a.timestep
            c3 = (induced_velocity(a.wing_aero_discretisation, c4) +
                induced_velocity(a.wake_discretisation, c4)) * a.timestep
            new_tevrs[i] = VortexRing(c1, c2, c3, c4)
        end
    else # We can use where the old rings were convected.
        for i = 1 : length(new_tevrs)
            c1 = coords(a.wing_aero_discretisation[i][end].geometry)[2]
            c4 = coords(a.wing_aero_discretisation[i][end].geometry)[3]
            c2 = coords(old_tevs[i].geometry)[1]
            c3 = coords(old_tevs[i].geometry)[4]
            new_tevrs[i] = VortexRing(c1, c2, c3, c4)
        end
    end
    a.last_shed_tevrs = new_tevrs
    # And turn the old ones into vortex particles to add to the wake.
    for ring in old_tevrs
        for fil in convert(Vector{StraightVortexFilament}, ring)
            append!(a.wake, get_children(to_particles(fil,
                a.shed_particle_radius, a.shed_particle_kernal)))
        end
    end
    return
end

#= CONTINIOUS MATHS ----------------------------------------------------------=#
function vorticity_densities(::Type{Lautat3D}, ref_vel::Float64,
    x_pos::S, coeffs::Vector{T}) where {S <: Real, T <: Real}

    @assert(length(coeffs > 0))
    @assert(-1 <= x_pos <= 1)

    theta = acos(-x_pos)
    mult = 2 * ref_vel
    ret = Vector{Float64}(undef, length(coeffs))
    ret[1] = mult * coeffs[i] * (1 + cos(theta)) / sin(theta)
    for i = 1 : length(coeffs)
        ret[i] = mult * coeffs[i] * sin((i-1)*theta)
    end
    return ret
end

function vorticity_density_derivatives(::Type{Lautat3D}, ref_vel::Float64,
    x_pos::S, coeffs::Vector{T}) where {S <: Real, T <: Real}

    reduced_fn = x->vorticity_densities(Lautat3D, ref_vel, x, coeffs)
    return ForwardDiff.derivative(reduced_fn, x_pos)
end

#= VECTOR-PROBLEM MAPPING ----------------------------------------------------=#
function fourier_coefficients_to_vector(a::Lautat3D)
    b = Vector{Float64}(undef, a.num_fourier_coeffs * length(a.strip_centres))
    @assert(length(a.strip_fourier_coeffs) == length(a.strip_centres),
     "length(strip_fourier_coeffs) != length(strip_centres).")
    for i = 1 : length(a.strip_centres)
        @assert(length(a.strip_fourier_coeffs[i]) == a.num_fourier_coeffs,
            string("Found incorrect number of fourier coefficents: ",
            "length(strip_fourier_coeffs[i]) == a.num_fourier_coeffs." ))
        b[(i - 1)*a.num_fourier_coeffs + 1 : i * a.num_fourier_coeffs] =
            a.strip_fourier_coeffs[i]
    end
    return b
end

function vector_to_fourier_coefficients!(
    a::Lautat3D, v::Vector{T}) where T <: Real

    @assert(length(v) == a.num_fourier_coeffs * length(a.strip_centres),
        string("Expected length of vector v to match n_fourier_coeffs *",
            "number of strips. length(v) was ", length(v), ", n_fourier_coeffs",
            " was ", a.num_fourier_coeffs, " and number of strips was ",
            length(a.strips_centres)))
    if any(isnan.(v)) @warn "NaN found in input vector!" end

    for i = 1 : length(a.strip_centres)
        a.strip_fourier_coeffs[i] =
            v[(i-1) * a.num_fourier_coeffs + 1 : i * a.num_fourier_coeffs]
    end
    return
end

function update_using_solution!(a::Lautat3D, solution::Vector{Float64},
    fourier_to_vort::Matrix{Float64})
    # The wing surface
    vorts = fourier_to_vort * solution
    update_using_vorticity_vector!(a.wing_aero_discretisation, vorts)
    # And the tev rings
    rings_per_strip = length(a.strip_discretisation)
    te_ring_indexes = collect(rings_per_strip: rings_per_strip :
                                    size(ring_vorticity_fourier_trans_mtrx, 1))
    update_using_vorticity_vector!(a.last_shed_tevrs, vorts[te_ring_indexes])
    return
end

#= UTILITY FUNCTIONS ---------------------------------------------------------=#
function get_low_and_high(a::Lautat3D, idx::Int, tipex::Bool=false)
    strip_pos = idx / (2 * a.strip_divisions) + 0.5
    strip_upper = Int64(ceil(strip_pos))
    strip_lower = Int64(floor(strip_pos))
    if tipex
        strip_upper = strip_upper > a.strip_divisions ?
            -1 : strip_upper
        strip_lower = strip_lower == 0 ? -1 : strip_lower
    else
        strip_upper = strip_upper > a.strip_divisions ?
            strip_lower : strip_upper
        strip_lower = strip_lower == 0 ? 1 : strip_lower
    end
    return (strip_lower, strip_upper, strip_pos % 1)
end

function foil_plane_normal(a::Lautat3D, strip_idx::Int, strip_yidx::Int)
    # Backward difference (appart from the first) to get foil tangent
    j = strip_yidx
    jp, cp = j > 1 ? (j-1, 1) : (j+1, -1)
    tangent = cp * unit(strip_eval_points[strip_idx][j] -
        a.strip_eval_points[strip_idx][jp])
    normal = unit(cross(tangent, a.strip_eval_normals[strip_idx][j]))
    return normal
end

function induced_normal_velocity_influence_vecmat(
    a::Lautat3D, inducing_body::Vorticity3D)
    # FIRST construct a big vector of the measurement points and normals...
    big_mes_pnt_vec, big_normal_vec = big_measurement_and_normal_vectors(a)
    # So that we can the map this to the influence function...
    point_influences = map(
        x->induced_velocity(inducing_body, x), big_mes_pnt_vec)
    # And the take only the influence normal to our surface.
    normal_influence = map(
        x->dot(x[1], x[2]), zip(point_influences, big_normal_vec))
    # Finally, we can turn our Vector{Vector} into a matrix.
    ret = zeros(length(normal_influence), length(normal_influence[1]))
    for row in 1 : size(ret, 1)
        b = normal_influence[row]
        println(b)
        println(ret[row, :])
        ret[row, :] = b
    end
    return ret
end

function big_measurement_and_normal_vectors(a::Lautat3D)
    @assert(length(a.strip_eval_points) == length(a.strip_centres),
        string("length(strip_eval_points) != length(strip_centres). ",
        "Perhaps strip_eval_points wasn't initalised? Try calling ",
        "generate_eval_points_data!(::Lautat3D) every time the ",
        "geometry changes. length(strip_eval_points) = ",
        length(a.strip_eval_points), ", lenght(strip_centres) = ",
        length(a.strip_centres)))
    stride = length(a.strip_eval_points[1])
    len = length(a.strip_eval_points) * length(a.strip_eval_points[1])
    big_mes_pnt_vec = Vector{Vector3D}(undef, len)
    big_normal_vec = Vector{Vector3D}(undef, len)
    for i = 1 : length(a.strip_eval_points)
        big_mes_pnt_vec[(i-1)*stride + 1 : i * stride] = a.strip_eval_points[i]
        big_normal_vec[(i-1)*stride + 1 : i * stride] = a.strip_eval_normals[i]
    end
    return big_mes_pnt_vec, big_normal_vec
end

#= MAPPING FOURIER TERMS TO RING STRENGTHS -----------------------------------=#
function ring_slant_corrections(a::Lautat3D)
    # The return type maps to the a.wing_aero_discretisation container.
    ret = Vector{Vector{Float64}}(undef, length(a.wing_aero_discretisation))

    for i = 1 : length(a.wing_aero_discretisation)
        ret[i] = zeros(length(a.strip_discretisation) - 1)
        #indexes of above and below strips.
        strip_upper, strip_lower, relpos = get_low_and_high(a, i)
        for j = 1 : length(a.strip_discretisation) - 1
            vring  = a.wing_aero_discretisation[i][j]
            # Foil plane normal y - plus/minus
            fpnyp = foil_plane_normal(a, strip_upper, j)
            fpnym = foil_plane_normal(a, strip_lower, j)
            # Generate a weighed of these normals according to the strip_pos
            fpnave = unit(relpos * unit(fpnyp) + (1 - relpos) * unit(fpnym))
            # Average LE and TE filament directions on ring:
            spanwise_dir = unit(unit(vring.c3 - vring.c2) +
                                                    unit(vring.c1 - vring.c4))
            ret[i][j] = 1 / dot(spanwise_dir, fpnave)
            if !isfinite(ret[i][j] ) @warn "Non-finite slant correction!" end
        end
    end
end


# The interaction matrix is controlled by a vector of ring vorticities.
# We need to generate the vector of ring vorticities from a vector of
# fourier terms. This is done via a transformation matrix.
function ring_vorticity_fourier_transformation_matrix(a::Lautat3D)
    # We have to contend with:
        # 1: Getting the value for the ring vorticity itself
        # 2: Correcting for the slant of the ring
        # 3: Interpolating between the strips on which we're doing the analysis.
    # Data structures should have a 1:1 index based mapping to the aero
    # discretisation.
    slant_corrections = ring_slant_corrections(a)

    centres_on_strips = a.strip_discretisation[1:end-1] +
        a.strip_discretisation[2:end]
    fourier_ones = ones(a.num_fourier_coeffs)
    # We get the strength due to each coefficient separately so Vector{Vector}
    generic_ring_strengths = map(
        x->vorticity_density_derivatives(Lautat3D, ref_vel, x, fourier_ones),
        centres_on_strips   )

    # We need an N_RINGS * N_FOURIER_HANDLES size matrix.
    n_per_strip = length(centres_on_strips)
    transform_matrix = zeros(n_per_strip * 2 * strip_divisions *
        length(a.strip_centres), length(a.strip_centres) * a.num_fourier_coeffs)

    for i = 1 : length(a.wing_aero_discretisation)
        strip_lower, strip_upper, relpos = get_low_and_high(a, i, true)
        for j = 1 : n_per_strip
            # The transform for a single ring is one row of the matrix.
            idx_ring = j + (i-1)*n_per_strip
            idxrange_fupper = collect(a.num_fourier_coeffs*(strip_upper-1)+1 :
                a.num_fourier_coeffs*strip_upper)
            idxrange_flower = collect(a.num_fourier_coeffs*(strip_lower-1)+1 :
                a.num_fourier_coeffs*strip_lower)
            if strip_upper != -1
            transform_matrix[idx_ring, idxrange_fupper] +=
                generic_ring_strengths[j] .* slant_corrections[i][j] .* relpos
            end
            if strip_lower != -1
            transform_matrix[idx_ring, idxrange_flower] +=
                generic_ring_strengths[j] .*slant_corrections[i][j] .*(1-relpos)
            end
        end
    end
    return transform_matrix
end

#= INFLUENCE MATRIX GENERATION -----------------------------------------------=#

# A matrix/vector of normal velocities induced on a surface must be converted
# to fourier coefficients via some kind of matrix. This will produce said
# matrix.
function induced_normal_velocity_to_fourier_transformation_matrix(
    a::Lautat3D, ref_vel :: Float64)

    eval_xs = (a.strip_discretisation[1:end-1]
        + a.strip_discretisation[2:end])/2
    eval_thetas = acos.(-eval_xs)
    eval_dthetas = acos.(-a.strip_discretisation[2:end]) .-
        acos.(-a.strip_discretisation[1:end-1])

    ret = zeros(length(a.strip_centres) * a.num_fourier_coeffs, length(eval_xs))
    for i = 1:size(ret, 1)
        j = i % a.num_fourier_coeffs
        c_terms = cos.(j .* eval_thetas)
        ret[i, :] = c_terms .* eval_dthetas .* (2   / (pi * ref_vel))
        if j == 0
            ret[i, :] = ret[i, :] .* -0.5
        end
    end
    return ret
end

# Compute the wing induced normal velocities matrix
function self_normal_velocity_influence_matrix(a::Lautat3D)
    @assert(length(a.wing_aero_discretisation) + 1 ==
        length(a.extended_strip_centers), string("Length of the ",
        "aero discretisation does not correspond to the strip discretisation.",
        " length(extended_strip_centres) - 1 != length(wing_aero_discretisat",
        "ion). They are ", length(a.extended_strip_centers) - 1, " and ",
        length(a.wing_aero_discretisation), " respectively."))
    ret =
        induced_normal_velocity_influence_vecmat(a, a.wing_aero_discretisation)
    return ret
end

function kinematics_normal_velocity_influence_vector(a::Lautat3D)
    point_influences = a.strip_eval_points_velocities
    # We can't use induced_normal_velocity_influence_vecmat(..) here,
    # but we can copy & paste most of it.
    big_mes_pnt_vec, big_normal_vec = big_measurement_and_normal_vectors(a)
    normal_influence = map(
        x->dot(x[1], x[2]), zip(point_influences, big_normal_vec))
    # Finally, we can turn our Vector{Vector} into a matrix.
    ret = zeros(length(normal_influence), length(normal_influence[1]))
    for row in 1 : size(ret, 1)
        ret[row, :] = normal_influence[i]
    end
    return ret
end

function psuedo_two_dimensional_induced_normal_velocity_matrix(a::Lautat3D)
    # This will be the psuedo 2D effect at each strip. Each strip only
    # affects itself here, so the matrix is block diagonal. The matrix
    # is also directly multiplied by the fourier term vector.
    cstride = length(a.strip_centres)
    ret = zeros(cstride * (length(a.strip_discretisation) - 1),
                a.num_fourier_coeffs * cstride)
    # We will position the 2D vortices where the rings LE and TE vortex fils
    # are.
    xps = a.strip_discretisation
    function v2D(xvort, xmes)  # 2D Vortex singularity. Katz pg60.
        dx = xmes .- xvort
        common = 1 / (2 * pi * (dx[1]^2 + dx[2]^2))
        return [dx[2] * common, -dx[1] * common]
    end
    for i = 1 : cstride
        vtx_pnts = a.strip_eval_points[i]
        vtx_nmls = a.strip_eval_normals[i]
        eval_pnts = (vtx_pnts[1: end-1] + vtx_pnts[2:end]) / 2
        eval_nmls = unit.((vtx_nmls[1: end-1] + vtx_nmls[2:end]) / 2)
        vtx_dxs = vtx_pnts[1: end-1] - vtx_pnts[2:end]
        # We need to make the problem 2D here. We establish a plane to work in
        # using the LE as the origin, and definitining a coordinate system
        # where x dir is based on TE - LE, out of plane is take as the normal
        # to the foil at the centre (we hope that the foil is straight...) and
        # generate y using the x-(out of plane) cross product.
        origin = eval_pts[1]
        x_dir = unit(eval_pts[end] - eval_pts[1])
        fp_nml = foil_plane_normal(a, i, floor(length(xps)/2))
        ydir = unit(cross(fp_nml, foil_dir))
        # The 2D coordinates of our vortex points:
        x2d_vort = zeros(2, length(vtx_pnts))
        for j = 1 : length(vtx_pnts)
            x2d_vort[:, j] = [  dot(foil_dir, vtx_pnts[j] - origin),
                                dot(ydir, vtx_pnts[j] - origin)     ]
        end
        # The 2D coordinates of our evaluation points:
        x2d_eval = zeros(2, length(eval_pnts))
        for j = 1 : length(eval_pnts)
            x2d_eval[:, j] = [  dot(foil_dir, eval_pnts[j] - origin),
                                dot(ydir, eval_pnts[j] - origin)     ]
        end
        # Now we can compute 2D influences:
        inf_mtrx = zeros(length(eval_pnts), length(vtx_pnts))
        for j = 1 : length(vtx_pnts)
            inf_mtrx[j, :] = map(
                x->v2D(vtx_pnts[j], x[1])' * x[2],
                zip(eval_pnts, eval_nmls))
        end
        # Convert this to vortex ring strength control
        vtx_ring_strength_transform = zeros(length(vtx_pnts), length(eval_pnts))
        vtx_ring_strength_transform[1:end-1,:] += LinearAlgebra.I
        vtx_ring_strength_transform[2:end,:] -= LinearAlgebra.I
        inf_mtrx *= vtx_ring_strength_transform
        # And now to fourier vector control.
        fourier_derivs = map(x->vorticity_density_derivatives(Lautat3D, ref_vel,
            x, ones(1:a.n_fourier_coeffs)),
            a.strip_discretisation[2:end] - a.strip_discretisation[1:end-1] )
        to_fourier_mtrx = zeros(length(eval_pnts), length(a.n_fourier_coeffs))
        for n = 1 : length(eval_pnts)
            to_fourier_mtrx[j, :] = fourier_derivs[j] .* vtx_dxs[j]
        end
        inf_matrix *= to_fourier_mtrx
        # We can leave the computing the fourier term effect to the transform
        # matrix we use elsewhere.
        ret[1 + (i-1) * a.num_fourier_coeffs : i * a.num_fourier_coeffs,
            1 + (i-1) * length(a.strip_discretisation) : i *
            length(strip_discretisation)] =
                inf_matrix
    end
    return inf_matrix
end

function tevr_fourier_induced_normal_velocities_matrix(
    a::Lautat3D,
    ring_vorticity_fourier_trans_mtrx::Matrix{Float64})
    # We want to know the trailing edge vortex ring strengths arising
    # from a given fourier vector.
    rings_per_strip = length(a.strip_discretisation)
    te_ring_indexes = collect(rings_per_strip: rings_per_strip :
                                    size(ring_vorticity_fourier_trans_mtrx, 1))
    te_ring_transforms = ring_vorticity_fourier_trans_mtrx[te_ring_indexes, :]
    # So that we can apply them to the incluence matrix we're about to make.
    ret = induced_normal_velocity_influence_vecmat(a, a.last_shed_tevrs)
    # And multiply through by a matrix to make sure the ring strenght is right.
    ret = ret * te_ring_transforms
    return ret
end

function solve_current_step(a::Lautat3D)
    integration_matrix =
        induced_normal_velocity_to_fourier_transformation_matrix(a, 1.0)
    # The "known" vector:
    wake_effect = induced_normal_velocity_influence_vecmat(a, a.wake)
    kinem_effect = kinematics_normal_velocity_influence_vector(a)
    known = integral_matrix * (wake_induced + kinematics_induced)

    # The "unknown" matrix
    fourier_to_ringstr = ring_vorticity_fourier_transformation_matrix(a)
    wing_effect = self_normal_velocity_influence_matrix(a)
    psuedo_2D_effect = psuedo_two_dimensional_induced_normal_velocity_matrix(a)
    tevr_effect = tevr_fourier_induced_normal_velocities_matrix(a,
        fourier_to_ringstr)
    unknown_matrix = integration_matrix * (LinearAlgebra.Indentity -
        integral_matrix *
        (wing_effect * fourier_to_ringstr + tevr_effect - psuedo_2D_effect))

    fourier_coeffs = unknown_matrix \ fourier_to_ringstr
    return fourier_coeffs
end

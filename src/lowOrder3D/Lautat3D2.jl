
mutable struct Lautat3D
    # Physical definitions:
    wing :: EquationSurf    # x chordwise (-1 = LE), y spanwise
    kinematics :: CoordinateTransform3D

    timestep :: Float64
    vortex_type :: Vortex3DRegularisationFunctions
    strip_centres :: Vector{Float64}
    strip_divisions :: Int64
    strip_discretisation :: Vector{Float64}
    num_fourier_coeffs :: Int64

    # Problem solving state:
    wake :: Vorticity3DSimpleCollector
    wing_geometry :: Vector{Vector{BilinearQuad}}
    wing_aero_discretisation :: Vorticity3DSimpleCollector
    # We have to have geometry each side of the strip. We have
    # length(strip_centres) * 2 * strip_divisions strips for this.
    extended_strip_centers :: Vector{Float64}
    strip_eval_points :: Vector{Vector{Vector3D}}
    strip_eval_normals :: Vector{Vector{Vector3D}}
    strip_eval_points_velocities :: Vector{Vector{Float64}}
    strip_fourier_coeffs :: Vector{Vector{Float64}}

    function Lautat3D()
        m_wing = EquationSurf(x->Vector3D(x[1], x[2], 0))
        m_kinematics = CoordinateTransform3D()
        timestep = NaN
        m_vortex_type = threed_planetary_kernels()
        m_strip_centres = Float64[]
        strip_divisions = 0
        strip_discretisation = []
        num_fourier_coeffs = 0
        wake = Vorticity3DSimpleCollector()
        wing_geometry = Vector{Vector{BilinearQuad}}()
        wing_aero_discretisation = Vorticity3DSimpleCollector()
        extended_strip_centers = Float64[]
        strip_eval_points = Vector{Vector{Vector3D}}()
        strip_eval_points_velocities = Vector{Vector{Vector3D}}()
        strip_eval_normals = Vector{Vector{Vector3D}}()
        strip_fourier_coeffs = Vector{Vector{Float64}}()
        new(m_wing, m_kinematics, timestep, m_vortex_type,
            m_strip_centres, strip_divisions, strip_discretisation,
            num_fourier_coeffs, wake, wing_geometry, wing_aero_discretisation,
            extended_strip_centers, strip_eval_points, strip_eval_normals,
            strip_fourier_coeffs)
    end
end

# Compute a.wing_geometry
function discretise_wing_to_geometry!(a::Lautat3D)
    @assert(length(a.strip_discretisation) > 20, string(
        "length(Lautat3D.strip_discretisation) probably needs to be more than",
        " 20.", " Current value is ", length(a.strip_discretisation), "."))
    @assert(extrema(a.strip_discretisation) == (-1, 1), string("All of the ",
        "Lautat3D.strip_discretisation values should be -1 <= val <= 1. ",
        "There should be points at -1 and 1. ",
        "strip_discretisation was ", a.strip_discretisation))
    @assert(issorted(a.strip_discretisation), "Expected "*
        "Lautat3D.strip_discretisation to be sorted in ascending order.")

    generate_extended_strip_centres!(a)
    a.wing_geometry = Vector{Vector{BilinearQuad}}(undef,
        2 * a.strip_divisions * length(a.strip_centres))
    for i = 1 : length(a.extended_strip_centers) - 1
        xgrid = a.strip_discretisation
        ygrid = [a.extended_strip_centers[i], a.extended_strip_centers[i+1]]
        a.wing_geometry[i] = discretise(a.wing, BilinearQuad, xgrid, ygrid)
    end
    return
end

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
    evl_x = (a.strip_discretisation[1:end-1] + a.strip_discretisation[2:end]) *
        0.5
    tmp_wing = EquationSurf(x->evaluate(a.kinematics, evaluate(a.wing, x)))
    for i = 1 : length(a.strip_centres)
        a.strip_eval_points = map(x->evaluate(tmp_wing, [x, a.strip_centres[i]]),
            evl_x)
        a.strip_eval_normals = map(x->unit(normal(tmp_wing,
            [x, a.strip_centres[i]])), evl_x)
        a.strip_eval_points_velocities = map(x->derivative(a.kinematics,
            evaluate(a.wing, [x, a.strip_centres])), evl_x)
    end
    return
end

function fourier_coefficients_to_vector(a::Lautat3D)
    b :: Vector{Float64}(undef, a.num_fourier_coeffs * length(a.strip_centres))
    @assert(length(a.strip_fourier_coeffs) == length(a.strip_centres),
     "length(strip_fourier_coeffs) != length(strip_centres).")
    for i = 1 : length(strip_centres)
        @assert(length(a.strip_fourier_coeffs[i]) == a.num_fourier_coeffs,
            string("Found incorrect number of fourier coefficents: ",
            "length(strip_fourier_coeffs[i]) == a.num_fourier_coeffs." ))
        b[(i - 1)*a.num_fourier_coeffs + 1 : i * a.num_fourier_coeffs] =
            a.strip_fourier_coeffs[i]
    end
    return b
end

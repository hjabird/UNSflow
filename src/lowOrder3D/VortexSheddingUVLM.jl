
mutable struct VortexSheddingUVLM
    # Inputs
    wing :: EquationSurf
    kinematics :: CoordinateTransform3D

    x_grid :: Vector{Float64}
    y_grid :: Vector{Float64}
    # Static state
    original_wing_discretisation :: Vector{BilinearQuad}
    # Time changing state:
    wing_geometry_discretisation :: Vector{BilinearQuad}
    wing_aero_discretisation :: Vorticity3DSimpleCollector
    wake :: Vorticity3DSimpleCollector
end

function discretise_wing!(a::VortexSheddingUVLM)
    @assert(length(x_grid) > 0)
    @assert(length(y_grid) > 0)

    geometry = discretise(a.wing, BilinearQuad, x_grid, y_grid)
    aero = Vorticity3DSimpleCollector
    for quad in geometry
        push!(aero, VortexRing(quad))
    end
    a.wing_geometry_discretisation = geometry
    a.wing_aero_discretisation = aero
    return
end

function wing_self_influence(a::VortexSheddingUVLM)
    centres = map(centre, a.wing_aero_discretisation)
    normals = map(normal, a.wing_geometry_discretisation)
    influence = map(
        x->vorticity_vector_velocity_influence(a.wing_aero_discretisation, x),
        centres)
    normal_influence = map(
        x->convert(Vector{Float}, x[1]) * x[2], zip(normals, incluence) )
    return normal_influence
end

function kinematics_influence(a::VortexSheddingUVLM)
    centres = map(centre, a.wing_aero_discretisation)
    normals = map(normal, a.wing_geometry_discretisation)
    influence = map(
        x-> dot(-derivative(a.kinem,x[1]), x[2]), zip(centres, normals))
    return normal_influence
end

function wake_influence(a::VortexSheddingUVLM)
    centres = map(centre, a.wing_aero_discretisation)
    normals = map(normal, a.wing_geometry_discretisation)
    ind_vels = map(x->induced_velocity(a.wake, x), centres)
    influence = map(x->dot(x[1], x[2]), zip(ind_vels, normals))
    return normal_influence
end

function solve_for_ring_strengths(a::VortexSheddingUVLM)
    mat = wing_self_influence(a)
    kinem_vect = kinematics_influence(a)
    wake_vect = wake_influence(a)
    new_vorticity = mat \ -(kinem_vect + wake_vect)
    update_using_vorticity_vector!(a.wing_aero_discretisation, new_vorticity)
    return
end

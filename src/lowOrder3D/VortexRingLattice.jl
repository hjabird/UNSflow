#===============================================================================
    VortexRingLattice.jl

    A vortex ring lattice object - reduces computational overhead compared
    to discrete vortex rings. This is NOT a Vorticity3DCollector since it 
    doesn't have a children field.

    Initial code: HJAB 2018

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

mutable struct VortexRingLattice <: Vorticity3D
    vorticity :: Matrix{Float64}
    geometry :: BilinearQuadSurf

    function VortexRingLattice(x_rings::T, y_rings::S) where {T<:Int, S<:Int}
        strs = zeros(x_rings, y_rings)
        geometry = BilinearQuadSurf(x_rings, y_rings)
        new(strs, geometry)
    end
    function VortexRingLattice(
        vorticity :: Matrix{Float64},
        geometry :: BilinearQuadSurf)

        @assert(size(vorticity) == size(geometry),
            string("Size of vorticity != Size of geometry.",
            " size(vorticity) = ", size(vorticity),
            " and size(geometry) = ", size(geometry), "."))
        if any(isnan.(vorticity)) || any(isinf.(vorticity))
            @warn "vorticity contained NaN or Inf value(s)."
        end
        new(vorticity, geometry)
    end
end

function VortexRingLattice(
    geometry :: BilinearQuadSurf)

    ni, nj = size(geometry)
    str = ones(ni, nj)
    return VortexRingLattice(str, geometry)
end

#= VORTICITY3D functions -----------------------------------------------------=#
function centre(a::VortexRingLattice)
    c = coords(a.geometry)
    sum = mapreduce(x->x, +, c)
    sum /= length(c)
    return sum
end

function effective_radius(a::VortexRingLattice)
    middle = centre(a)
    rad = mapreduce(x->abs(middle - x), max, coords(a.geometry))
    return rad
end

function vorticity(a::VortexRingLattice)
    return Vector3D(0,0,0)
end

function induced_velocity(a::VortexRingLattice, measurement_loc::Vector3D)
    vel = Vector3D(0,0,0)
    if length(a.vorticity) == 0
        return vel
    end
    id, jd = size(a.vorticity)
    c = a.geometry.coordinates
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start = c[i, j]
            stop = c[i, j + 1]
            str = a.vorticity[i - 1, j] - a.vorticity[i, j]
            vel += induced_velocity(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.vorticity[i, j] - a.vorticity[i, j - 1]
            vel += induced_velocity(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # external horizontal filaments
    imax, jmax = size(c)
    vel = mapreduce(
        i->induced_velocity(StraightVortexFilament, c[i, 1],
            c[i+1, 1], a.vorticity[i, 1], measurement_loc),
        +, 1 : id; init=vel)
    vel = mapreduce(
        i->induced_velocity(StraightVortexFilament, c[i, jmax],
            c[i+1, jmax], -a.vorticity[i, jd], measurement_loc),
        +, 1 : id; init=vel)
    vel = mapreduce(
        j->induced_velocity(StraightVortexFilament, c[1, j],
            c[1, j+1], -a.vorticity[1, j], measurement_loc),
        +, 1 : jd; init=vel)
    vel = mapreduce(
        j->induced_velocity(StraightVortexFilament, c[imax, j],
            c[imax, j+1], a.vorticity[id, j], measurement_loc),
        +, 1 : jd; init=vel)
    return vel
end

function induced_velocity_curl(a::VortexRingLattice)
    curl = zeros(3,3)
    if length(a.vorticity) == 0
        return curl
    end
    id, jd = size(a.vorticity)
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start =c[i, j]
            stop = c[i, j + 1]
            str = a.vorticity[i - 1, j] - a.vorticity[i, j]
            curl += induced_velocity_curl(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.vorticity[i, j] - a.vorticity[i, j - 1]
            curl += induced_velocity_curl(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # external horizontal filaments
    imax, jmax = size(c)
    curl = mapreduce(
        i->induced_velocity_curl(StraightVortexFilament, c[i, 1],
        c[i+1, 1], a.vorticity[i, 1], measurement_loc),
        1 : id; init=curl)
    curl = mapreduce(
        i->induced_velocity_curl(StraightVortexFilament, c[i, jmax],
            c[i+1, jmax], -a.vorticity[i, jd], measurement_loc),
        1 : id; init=curl)
    curl = mapreduce(
        j->induced_velocity_curl(StraightVortexFilament, c[1, j],
        c[1, j+1], -a.vorticity[1, j], measurement_loc),
        1 : jd; init=curl)
    curl = mapreduce(
        j->induced_velocity_curl(StraightVortexFilament, c[imax, j],
        c[imax, j+1], a.vorticity[id, j], measurement_loc),
         1 : jd; init=curl)
    return curl
end

function steady_forces(a::VortexRingLattice,
    vel_fn::Function, 
    density::Real=1;
    imax_filament_strs::Vector{T}=Vector{Float64}(),   # Bottom edge
    jmax_filament_strs::Vector{T}=Vector{Float64}(),   # Right edge
    i1_filament_strs::Vector{T}=Vector{Float64}()  ,   # Top edge
    j1_filament_strs::Vector{T}=Vector{Float64}()  ) where T <: Real # Left edge
    # INPUT CHECKS
    @assert(hasmethod(vel_fn, Tuple{Vector3D}),
        "Expected induced_vel_fn passed to steady force function to take a "*
        "single arguement of type UNSflow.Vector3D.")
    function asserttmplt(argname, arg, dim)
        @assert((length(arg) == 0) || (length(arg) ==
            size(a.vorticity, dim)), string("Incorrect length vector ",
            "passed through keyword argument ", argname, ". Length was ",
            length(arg), ". Should have been ", 
            size(a.vorticity, dim), "."))
        if any(isnan.(arg)) || any(isinf.(arg))
            @warn string("Vector passed through keyword arguement ",
                argname, " has NaN or Inf compents.")
        end
    end
    asserttmplt("i1_filament_strs", i1_filament_strs, 2)
    asserttmplt("imax_filament_strs", imax_filament_strs, 2)
    asserttmplt("j1_filament_strs", j1_filament_strs, 1)
    asserttmplt("jmax_filament_strs", jmax_filament_strs, 1)
    replacefn = (x, dim)-> length(x) > 0 ? x : zeros(size(a.vorticity, dim))
    i1_filament_strs = replacefn(i1_filament_strs, 2)
    imax_filament_strs = replacefn(imax_filament_strs, 2)
    j1_filament_strs = replacefn(j1_filament_strs, 1)
    jmax_filament_strs = replacefn(jmax_filament_strs, 1)

    # AND the actual maths bit:
    id, jd = size(a.vorticity)
    c = a.geometry.coordinates
    forces = map(x->Vector3D(0,0,0), [0 for i = 1 : id, j = 1 : jd])    
    force_fn = (sta, sto, str)->steady_force(StraightVortexFilament, sta, sto, 
                                                        str, vel_fn, density)
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start = c[i, j]
            stop = c[i, j + 1]
            str = a.vorticity[i - 1, j] - a.vorticity[i, j]
            force = force_fn(start, stop, str)
            forces[i-1, j] += force / 2
            forces[i, j] += force / 2
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.vorticity[i, j] - a.vorticity[i, j - 1]
            force = force_fn(start, stop, str)
            forces[i, j-1] += force / 2
            forces[i, j] += force / 2
        end
    end
    # Now the external filaments.    
    imax, jmax = size(c)
    forces[id, :] += map(
        j->force_fn(c[imax, j], c[imax, j+1], 
            a.vorticity[id, j] + imax_filament_strs[j]), 
        collect(1:jmax-1))
    forces[:, jd] += map(
        i->force_fn(c[i+1, jmax], c[i, jmax], 
            a.vorticity[i, jd] + jmax_filament_strs[id - i + 1]), 
        collect(1:imax-1))
    forces[1, :] += map(
        j->force_fn(c[1, j+1], c[1, j], 
            a.vorticity[1, j] + i1_filament_strs[jd - j + 1]), 
        collect(1:jmax-1))
    forces[:, 1] += map(
        i->force_fn(c[i, 1], c[i + 1, 1], 
            a.vorticity[i, 1] + j1_filament_strs[i]), 
        collect(1:imax-1))
    return forces
end

function unsteady_forces(
    a::VortexRingLattice, 
    old_vort_vect::Vector{T},
    dt::Real) where T <: Real

    # Unsteady Bernoulli:
    # (P / rho) = Q^2/2 - vref^2/2 + dpotential/dt
    # (P / rho) = steady_result + dpotential/dt
    # So we want the difference in potential over the top/bottom of the ring.
    # This is approximately the rate of change of ring vorticity
    # See Katz 13.149
    @assert(vorticity_vector_length(a) == length(old_vort_vect),
        string("Expected old and current vorticity vector lengths to match.",
        " Old vector length is ", length(old_vort_vect), " and current ",
        "vector length is ", vorticity_vector_length(a)))
    @assert(dt != 0, string("dt cannot be equal to zero. Given value was ",
                            dt))
    if any(isnan.(old_vort_vect)) || any(isinf.(old_vort_vect))
        @warn "Inf or NaN values in old vorticity vector"
    end
    current_vort_vector = vorticity_vector(a)
    if any(isnan.(current_vort_vector)) || any(isinf.(current_vort_vector))
        @warn "Inf or NaN values in current vorticity vector"
    end
    dvort = (vorticity_vector(a) - old_vort_vect) / dt
    m_unsteady_pressures = dvort
    m_areas = vec(areas(a.geometry))
    m_normals = vec(normals(a.geometry))
    m_forces = map(x->x[1]*unit(x[2])*x[3], 
        zip(m_unsteady_pressures, m_normals, m_areas))
    m_forces = reshape(m_forces, size(a.vorticity))
    return m_forces
end

function steady_pressures(a::VortexRingLattice,
    vel_fn::Function, 
    density::Real=1;
    imax_filament_strs::Vector{T}=Vector{Float64}(),   # Bottom edge
    jmax_filament_strs::Vector{T}=Vector{Float64}(),   # Right edge
    i1_filament_strs::Vector{T}=Vector{Float64}()  ,   # Top edge
    j1_filament_strs::Vector{T}=Vector{Float64}()  ) where T <: Real # Left edge

    replacefn = (x, dim)-> length(x) > 0 ? x : zeros(size(a.vorticity, dim))
    i1_filament_strs = replacefn(i1_filament_strs, 2)
    imax_filament_strs = replacefn(imax_filament_strs, 2)
    j1_filament_strs = replacefn(j1_filament_strs, 1)
    jmax_filament_strs = replacefn(jmax_filament_strs, 1)

    force = steady_forces(a, vel_fn, density;
        imax_filament_strs=imax_filament_strs,
        jmax_filament_strs=jmax_filament_strs,
        i1_filament_strs=i1_filament_strs,
        j1_filament_strs=j1_filament_strs)
    mareas = areas(a.geometry)
    mnormals = normals(a.geometry)
    normal_forces = map(x->dot(x[1], x[2]), zip(mnormals, force))
    press = normal_forces ./ mareas
    return press
end

function unsteady_pressures(a::VortexRingLattice,
    old_vort_vect::Vector{T}, dt::Real) where T <: Real

    forces = unsteady_forces(a, old_vort_vect, dt)
    m_areas = areas(a.geometry)
    m_normals = normals(a.geometry)
    m_press = map(x->dot(x[1], unit(x[2]))/x[3], 
        zip(forces, m_normals, m_areas))
    return m_press
end

function steady_loads(a::VortexRingLattice,
    velocity_fn::Function, 
    density::Real=1;
    measurement_centre::Vector3D=Vector3D(0,0,0),
    imax_filament_strs::Vector{T}=Vector{Float64}(),   # Bottom edge
    jmax_filament_strs::Vector{T}=Vector{Float64}(),   # Right edge
    i1_filament_strs::Vector{T}=Vector{Float64}()  ,   # Top edge
    j1_filament_strs::Vector{T}=Vector{Float64}()  ) where T <: Real # Left edge

    replacefn = (x, dim)-> length(x) > 0 ? x : zeros(size(a.vorticity, dim))
    i1_filament_strs = replacefn(i1_filament_strs, 2)
    imax_filament_strs = replacefn(imax_filament_strs, 2)
    j1_filament_strs = replacefn(j1_filament_strs, 1)
    jmax_filament_strs = replacefn(jmax_filament_strs, 1)

    mforces = steady_forces(a, velocity_fn, density;
        imax_filament_strs=imax_filament_strs,
        jmax_filament_strs=jmax_filament_strs,
        i1_filament_strs=i1_filament_strs,  
        j1_filament_strs=j1_filament_strs)
    mcents = centres(a.geometry)
    total_forces = sum(mforces)
    moment_contrib = map(
        x->cross(x[1]-measurement_centre, x[2]),
        zip(mcents, mforces))
    total_moment = sum(moment_contrib)
    # Calculate pressures... meh...
    mareas = areas(a.geometry)
    mnormals = normals(a.geometry)
    normal_forces = map(x->dot(x[1], x[2]), zip(mnormals, mforces))
    all_pressures = normal_forces ./ mareas
    return total_forces, total_moment, all_pressures
end

function unsteady_loads(a::VortexRingLattice,
    old_vort_vect::Vector{T}, dt::Real;
    measurement_centre::Vector3D=Vector3D(0,0,0)) where T <: Real

    mpress = unsteady_pressures(a, old_vort_vect, dt)
    mforces = unsteady_forces(a, old_vort_vect, dt)
    mareas = areas(a.geometry)
    mcents = centres(a.geometry)
    total_forces = sum(mforces)
    moment_contrib = map(
        x->cross(x[1]-measurement_centre, x[2]),
        zip(mcents, mforces))
    total_moment = sum(moment_contrib)
    return total_forces, total_moment, mpress
end

function state_vector_length(a::VortexRingLattice)
    return 3 * length(a.geometry.coordinates)
end

function state_vector(a::VortexRingLattice)
    v = zeros(state_vector_length(a))
    for i = 1 : length(a.geometry.coordinates)
        v[(i-1) * 3 + 1 : i * 3] = convert(Vector{Float64}, 
            a.geometry.coordinates[i])
    end
    return v
end

function update_using_state_vector!(
    this::VortexRingLattice,
    state_vect::Vector{Float64})

    n = length(coords(this.geometry))
    @assert(length(state_vect) == n * 3, string(
        "State vector was incorrect length. Length was ", length(state_vect),
        " but should have been ", n * 3, "."))
    if any(isnan.(state_vect)) || any(isinf.(state_vect))
        @warn "Infinite or NaN values in state vector."
    end
    for i = 1 :  n
        this.geometry.coordinates[i] = 
            convert(Vector3D, state_vect[i*3 - 2: i*3])
    end
    return
end

function state_time_derivative(
    this::VortexRingLattice,
    inducing_bodies::Vorticity3D)

    svtd = Vector{Float64}(undef, 3 * length(this.geometry.coordinates))
    c = coords(this.geometry)
    for i = 1 : length(c)
        svtd[i*3 - 2 : i*3] = 
            convert(Vector{Float64}, induced_velocity(inducing_bodies, c[i]))
    end
    return svtd
end

function vorticity_vector_length(this::VortexRingLattice)
    return length(this.vorticity)
end

function vorticity_vector(this::VortexRingLattice)
    return deepcopy(vec(this.vorticity))
    # Otherwise updating this.vorticity will pdate the output vorticity_vector
    # which is no use if you were trying to store the old one for example.
end

function update_using_vorticity_vector!(
    this::VortexRingLattice, vort_vect::Vector{Float64})

    @assert(length(vort_vect) == vorticity_vector_length(this),
        string("Incorrect vorticity_vector length for VortexRingLattice.",
        " Expected a length of ", vorticity_vector_length(this), " but ",
        " was given a vector of length ", length(vort_vect)))
    for i = 1 : length(vort_vect)
        this.vorticity[i] = vort_vect[i]
    end
    return
end

function vorticity_vector_velocity_influence(
    this :: VortexRingLattice,
    mes_pnt :: Vector3D)

    ret = zeros(3, vorticity_vector_length(this))
    indexes = [(i, j) for 
        i = 1 : size(this.vorticity, 1), j = 1 : size(this.vorticity, 2)]
    indexes = vec(indexes)
    for i = 1 : length(indexes)
        idx = indexes[i]
        quad = convert(BilinearQuad, this.geometry, idx[1], idx[2])
        ring = VortexRing(quad, 1.)
        influence = induced_velocity(ring, mes_pnt)
        ret[:, i] = convert(Vector{Float64}, influence)
    end
    return ret
end

function extract_vorticity_vector_indexes(
    a::VortexRingLattice,
    i_ind :: Vector{Int},
    j_ind :: Vector{Int})
    
    ni, nj = size(a.vorticity)
    for i = 1 : length(i_ind)
        @assert(0 < i_ind[i] <= ni, string("Tried to extract invalid i index",
        " from VortexRingLattice. Desired value ", i_ind[i], " at i_ind[",
        i, "] exceeded max index of ", ni, "."))
    end
    for j = 1 : length(j_ind)
        @assert(0 < j_ind[j] <= nj, string("Tried to extract invalid j index",
        " from VortexRingLattice. Desired value ", j_ind[j], " at j_ind[",
        j, "] exceeded max index of ", nj, "."))
    end
    lookup = [ni * (j-1) + i for i = 1 : ni, j = 1 : nj]
    ret = map(x->lookup[x[1], x[2]], [(i, j) for i in i_ind, j in j_ind])
    return ret
end

function add_row!(
    a::VortexRingLattice, 
    points::Vector{Vector3D},
    strengths::Vector{T}=zeros(1, size(a.coords, 2) - 1)) where T <: Real
    coords = a.geometry.coordinates
    @assert(size(coords, 2) == length(points), string(
        "Dimension mismatch: size(VortRingLattice.geometry.coords, 2) must ",
        "match length(points). length(points) = ", length(points), " and ",
        " size(VortRingLattice.geometry.coords, 2) = ", size(coords, 2)))
    a.geometry.coordinates = vcat(coords, reshape(points, 1, :))
    a.vorticity = vcat(a.vorticity, reshape(strengths, 1, :))
    return
end

#= CONVERSION FUNCTIONS ------------------------------------------------------=#
function Base.convert(::Type{Matrix{VortexRing}}, a::VortexRingLattice)
    ret = Matrix{VortexRing}(undef, size(a.vorticity))
    g = a.geometry.coordinates
    s = a.vorticity
    ret = map(
        x->VortexRing(g[x[1]+1, x[2]], g[x[1]+1, x[2]+1], g[x[1], x[2]+1],
            g[x[1], x[2]], s[x[1], x[2]] ),
        [(i, j) for i = 1 : size(s,1), j = 1 : size(s, 2)]
        )
    return ret
end

function Base.convert(::Type{Vector{VortexRing}}, a::VortexRingLattice)
    return vec(Base.convert(Matrix{VortexRing}, a))
end

function Base.convert(
    ::Type{Vector{StraightVortexFilament}}, 
    a::VortexRingLattice)

    id, jd = size(a.vorticity)
    c = a.geometry.coordinates
    fils = Vector{StraightVortexFilament}(undef, id * jd * 2 + id + jd)
    offset = 1
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start = c[i, j]
            stop = c[i, j + 1]
            str = a.vorticity[i - 1, j] - a.vorticity[i, j]
            fils[offset] = StraightVortexFilament(start, stop, str)
            offset += 1
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.vorticity[i, j] - a.vorticity[i, j - 1]
            fils[offset] = StraightVortexFilament(start, stop, str)
            offset += 1
        end
    end
    # external horizontal filaments
    imax, jmax = size(c)
    fils[offset : offset + id - 1] = map(
        i->StraightVortexFilament(c[i, 1], c[i+1, 1], a.vorticity[i, 1]),1 : id)
    offset += id
    fils[offset : offset + id - 1] = map(
        i->StraightVortexFilament(c[i, jmax], c[i+1, jmax],-a.vorticity[i, jd]),
        1 : id)
    offset += id
    fils[offset : offset + jd - 1] = map(
        j->StraightVortexFilament(c[1, j], c[1, j+1], -a.vorticity[1, j]),
        1 : jd)
    offset += jd
    fils[offset : offset + jd - 1] = map(
        j->StraightVortexFilament(c[imax, j], c[imax, j+1], a.vorticity[id, j]),
        1 : jd)
    return fils
end

function Base.push!(a::UnstructuredMesh, b::VortexRingLattice, 
    controldict=Dict{String, Any}())

    if haskey(controldict, "VortexRingLatticeAsRings") &&
        controldict["VortexRingLatticeAsRings"] == true
        rings = convert(Vector{VortexRing}, b)
        for ring in rings
            push!(a, ring, controldict)
        end
    elseif haskey(controldict, "VortexRingLatticeAsFilaments") &&
        controldict["VortexRingLatticeAsFilaments"] == true

        fils = convert(Vector{StraightVortexFilament}, b)
        for fil in fils
            push!(a, fil, controldict)
        end
    else # DEFAULT
        quads = convert(Matrix{BilinearQuad}, b.geometry)
        idxs = [(i, j) for i = 1 : size(quads, 1), j = 1 : size(quads, 2)]
        for idx in idxs
            cidx, pidx = add_cell!(a, quads[idx[1], idx[2]])
            add_celldata!(a, cidx, "Filament_vorticity", 
                b.vorticity[idx[1], idx[2]])
        end
    end
    return
end

function add_celldata!(a::MeshDataLinker, b::VortexRingLattice, 
    dataname::String, data::Union{Matrix{Float64}, Matrix{Vector3D}})

    @assert(size(data) == size(b.geometry), 
        string("size(data) must equal the size(geometry) of the vortex",
            " lattice. size(data) = ", size(data), " and size(geometry) = ",
            size(b.geometry), "."))
    add_celldata!(a, b, dataname, vec(data))
    return
end

function add_celldata!(a::MeshDataLinker, b::VortexRingLattice, 
    dataname::String, data::Union{Vector{Float64}, Vector{Vector3D}})

    @assert(length(data) == length(b.geometry), 
        string("length(data) must equal the length(geometry) of the vortex",
            " lattice. length(data) = ",length(data)," and length(geometry) = ",
            length(b.geometry), "."))
    geom1 = convert(Vector{BilinearQuad}, b.geometry)
    geom2 = vec(map(x->x.geometry, convert(Vector{VortexRing}, b)))
    for i = 1 : length(geom1)
        add_celldata!(a, dataname, geom1[i], data[i])
        add_celldata!(a, dataname, geom2[i], data[i])
    end
    return
end

function add_pointdata!(a::MeshDataLinker, b::VortexRingLattice, 
    dataname::String, data::Union{Matrix{Float64}, Matrix{Vector3D}})

    @assert(size(data) == size(b.geometry.coordinates), 
        string("size(data) must equal the size(geometry.coordinates) of the",
            " vortex lattice. size(data) = ", size(data), 
            " and size(geometry) = ", size(b.geometry.coordinates), "."))
    add_pointdata!(a, b, dataname, vec(data))
end

function add_pointdata!(a::MeshDataLinker, b::VortexRingLattice, 
    dataname::String, data::Union{Vector{Float64}, Vector{Vector3D}})

    points = vec(b.geometry.coordinates)
    @assert(length(data) == length(b.coordinates), 
        string("The length of the data vector should was ", length(data),
        " which did not match length(VortexRingLattice.coordinates) =  ",
        length(b.geometry.coordinates)))
    for v in zip(data, points)
        add_pointdata!(a, dataname, v[2], v[1])
    end
    return
end

function Base.print(a::VortexRingLattice)
    print(typeof(a), " with size ", size(a.geometry), ". ")
end

function Base.println(a::VortexRingLattice)
    print(typeof(a), " with size ", size(a.geometry), ".\n ")
end

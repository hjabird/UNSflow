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
    strengths :: Matrix{Float64}
    geometry :: BilinearQuadSurf

    function VortexRingLattice(x_rings::Float64, y_rings::Float64)
        strs = zeros(x_rings, y_rings)
        geometry = BilinearQuadSurf(x_rings, y_rings)
        new(strs, geometry)
    end
    function VortexRingLattice(
        strengths :: Matrix{Float64},
        geometry :: BilinearQuadSurf)

        @assert(size(strengths) == size(geometry),
            string("Size of strengths != Size of geometry.",
            " size(strengths) = ", size(strengths),
            " and size(geometry) = ", size(geometry), "."))
        if any(isnan.(strengths)) || any(isinf.(strengths))
            @warn "strengths contained NaN or Inf value(s)."
        end
        new(strengths, geometry)
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
    id, jd = size(a.strengths)
    c = a.geometry.coordinates
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start = c[i, j]
            stop = c[i, j + 1]
            str = a.strengths[i - 1, j] - a.strengths[i, j]
            vel += induced_velocity(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.strengths[i, j] - a.strengths[i, j - 1]
            vel += induced_velocity(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # external horizontal filaments
    imax, jmax = size(c)
    vel = mapreduce(
        i->induced_velocity(StraightVortexFilament, c[i, 1],
            c[i+1, 1], a.strengths[i, 1], measurement_loc),
        +, 1 : id; init=vel)
    vel = mapreduce(
        i->induced_velocity(StraightVortexFilament, c[i, jmax],
            c[i+1, jmax], -a.strengths[i, jd], measurement_loc),
        +, 1 : id; init=vel)
    vel = mapreduce(
        j->induced_velocity(StraightVortexFilament, c[1, j],
            c[1, j+1], -a.strengths[1, j], measurement_loc),
        +, 1 : jd; init=vel)
    vel = mapreduce(
        j->induced_velocity(StraightVortexFilament, c[imax, j],
            c[imax, j+1], a.strengths[id, j], measurement_loc),
        +, 1 : jd; init=vel)
    return vel
end

function induced_velocity_curl(a::VortexRingLattice)
    curl = zeros(3,3)
    id, jd = size(a.strengths)
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start =c[i, j]
            stop = c[i, j + 1]
            str = a.strengths[i - 1, j] - a.strengths[i, j]
            curl += induced_velocity_curl(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.strengths[i, j] - a.strengths[i, j - 1]
            curl += induced_velocity_curl(StraightVortexFilament, 
                start, stop, str, measurement_loc)
        end
    end
    # external horizontal filaments
    imax, jmax = size(c)
    curl = mapreduce(
        i->induced_velocity_curl(StraightVortexFilament, c[i, 1],
        c[i+1, 1], a.strengths[i, 1], measurement_loc),
        1 : id; init=curl)
    curl = mapreduce(
        i->induced_velocity_curl(StraightVortexFilament, c[i, jmax],
            c[i+1, jmax], -a.strengths[i, jd], measurement_loc),
        1 : id; init=curl)
    curl = mapreduce(
        j->induced_velocity_curl(StraightVortexFilament, c[1, j],
        c[1, j+1], -a.strengths[1, j], measurement_loc),
        1 : jd; init=curl)
     curl = mapreduce(
        j->induced_velocity_curl(StraightVortexFilament, c[imax, j],
        c[imax, j+1], a.strengths[id, j], measurement_loc),
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
<<<<<<< HEAD
            size(a.strengths, dim)), string("Incorrect length vector ",
            "passed through keyword argument ", argname, ". Length was ",
            length(arg), ". Should have been ", 
            size(a.strengths, dim), "."))
=======
            size(a.coordinates, dim) - 1), string("Incorrect length vector ",
            "passed through keyword argument ", argname, ". Length was ",
            length(arg), ". Should have been ", 
            size(a.coordinates, dim) - 1, "."))
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
        if any(isnan.(arg)) || any(isinf.(arg))
            @warn string("Vector passed through keyword arguement ",
                argname, " has NaN or Inf compents.")
        end
    end
<<<<<<< HEAD
    asserttmplt("i1_filament_strs", i1_filament_strs, 2)
    asserttmplt("imax_filament_strs", imax_filament_strs, 2)
    asserttmplt("j1_filament_strs", j1_filament_strs, 1)
    asserttmplt("jmax_filament_strs", jmax_filament_strs, 1)
    replacefn = (x, dim)-> length(x) > 0 ? x : zeros(size(a.strengths, dim))
    i1_filament_strs = replacefn(i1_filament_strs, 2)
    imax_filament_strs = replacefn(imax_filament_strs, 2)
    j1_filament_strs = replacefn(j1_filament_strs, 1)
    jmax_filament_strs = replacefn(jmax_filament_strs, 1)
=======
    asserttmplt("i1_filament_strs", i1_filament_strs, 1)
    asserttmplt("imax_filament_strs", imax_filament_strs, 1)
    asserttmplt("j1_filament_strs", j1_filament_strs, 2)
    asserttmplt("jmax_filament_strs", jmax_filament_strs, 2)
    replacefn = (x, dim)-> length(x) > 0 ? x : size(a.strengths, dim)
    i1_filament_strs = replacefn(i1_filament_strs, 1)
    imax_filament_strs = replacefn(imax_filament_strs, 1)
    j1_filament_strs = replacefn(j1_filament_strs, 2)
    jmax_filament_strs = replacefn(jmax_filament_strs, 2)
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e

    # AND the actual maths bit:
    id, jd = size(a.strengths)
    c = a.geometry.coordinates
<<<<<<< HEAD
    forces = map(x->Vector3D(0,0,0), [0 for i = 1 : id, j = 1 : jd])    
=======
    forces = map(Vector3D(0,0,0), [0 for i = 1 : id, j = 1 : jd])    
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
    force_fn = (sta, sto, str)->steady_force(StraightVortexFilament, sta, sto, 
                                                        str, vel_fn, density)
    # internal horizontal filaments
    for i = 2 : id
        for j = 1 : jd
            start = c[i, j]
            stop = c[i, j + 1]
            str = a.strengths[i - 1, j] - a.strengths[i, j]
<<<<<<< HEAD
            force = force_fn(start, stop, str)
=======
            force = force_fn(start, strop, str)
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
            forces[i-1, j] += force / 2
            forces[i, j] += force / 2
        end
    end
    # internal vertical filaments
    for i = 1 : id
        for j = 2 : jd
            start = c[i, j]
            stop = c[i + 1, j]
            str = a.strengths[i, j] - a.strengths[i, j - 1]
<<<<<<< HEAD
            force = force_fn(start, stop, str)
=======
            force = force_fn(start, strop, str)
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
            forces[i, j-1] += force / 2
            forces[i, j] += force / 2
        end
    end
    # Now the external filaments.    
    imax, jmax = size(c)
    forces[id, :] += map(
        j->force_fn(c[imax, j], c[imax, j+1], 
            a.strengths[id, j] + imax_filament_strs[j]), 
        collect(1:jmax-1))
    forces[:, jd] += map(
        i->force_fn(c[i+1, jmax], c[i, jmax], 
<<<<<<< HEAD
            a.strengths[i, jd] + jmax_filament_strs[id - i + 1]), 
        collect(1:imax-1))
    forces[1, :] += map(
        j->force_fn(c[1, j+1], c[1, j], 
            a.strengths[1, j] + i1_filament_strs[jd - j + 1]), 
        collect(1:jmax-1))
    forces[:, 1] += map(
        i->force_fn(c[i, 1], c[i + 1, 1], 
            a.strengths[i, 1] + j1_filament_strs[i]), 
=======
            a.strengths[i, jd] + imax_filament_strs[id - i + 1]), 
        collect(1:imax-1))
    forces[1, :] += map(
        j->force_fn(c[1, j+1], c[1, j], 
            a.strengths[1, j] + imax_filament_strs[ij - j + 1]), 
        collect(1:jmax-1))
    forces[:, 1] += map(
        i->force_fn(c[i, 1], c[i + 1, 1], 
            a.strengths[i, 1] + imax_filament_strs[i]), 
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
        collect(1:imax-1))
    return forces
end

function steady_pressures(a::VortexRingLattice,
    vel_fn::Function, 
    density::Real=1;
    imax_filament_strs::Vector{T}=Vector{Float64}(),   # Bottom edge
    jmax_filament_strs::Vector{T}=Vector{Float64}(),   # Right edge
    i1_filament_strs::Vector{T}=Vector{Float64}()  ,   # Top edge
    j1_filament_strs::Vector{T}=Vector{Float64}()  ) where T <: Real # Left edge

<<<<<<< HEAD
    replacefn = (x, dim)-> length(x) > 0 ? x : zeros(size(a.strengths, dim))
    i1_filament_strs = replacefn(i1_filament_strs, 2)
    imax_filament_strs = replacefn(imax_filament_strs, 2)
    j1_filament_strs = replacefn(j1_filament_strs, 1)
    jmax_filament_strs = replacefn(jmax_filament_strs, 1)

    force = steady_forces(a, vel_fn, density;
=======
    replacefn = (x, dim)-> length(x) > 0 ? x : size(a.strengths, dim)
    i1_filament_strs = replacefn(i1_filament_strs, 1)
    imax_filament_strs = replacefn(imax_filament_strs, 1)
    j1_filament_strs = replacefn(j1_filament_strs, 2)
    jmax_filament_strs = replacefn(jmax_filament_strs, 2)

    force = steady_forces(a, vel_fn, density,
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
        imax_filament_strs=imax_filament_strs,
        jmax_filament_strs=jmax_filament_strs,
        i1_filament_strs=i1_filament_strs,
        j1_filament_strs=j1_filament_strs)
    mareas = areas(a.geometry)
<<<<<<< HEAD
    mnormals = normals(a.geometry)
    normal_forces = map(x->dot(x[1], x[2]), zip(mnormals, force))
    press = normal_forces ./ mareas
=======
    press = forces ./ mareas
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
    return press
end

function steady_loads(a::VortexRingLattice,
    velocity_fn::Function, 
    density::Real=1;
    measurement_centre::Vector3D=Vector3D(0,0,0),
    imax_filament_strs::Vector{T}=Vector{Float64}(),   # Bottom edge
    jmax_filament_strs::Vector{T}=Vector{Float64}(),   # Right edge
    i1_filament_strs::Vector{T}=Vector{Float64}()  ,   # Top edge
    j1_filament_strs::Vector{T}=Vector{Float64}()  ) where T <: Real # Left edge

<<<<<<< HEAD
    replacefn = (x, dim)-> length(x) > 0 ? x : zeros(size(a.strengths, dim))
    i1_filament_strs = replacefn(i1_filament_strs, 2)
    imax_filament_strs = replacefn(imax_filament_strs, 2)
    j1_filament_strs = replacefn(j1_filament_strs, 1)
    jmax_filament_strs = replacefn(jmax_filament_strs, 1)

    mforces = steady_forces(a, velocity_fn, density;
=======
    replacefn = (x, dim)-> length(x) > 0 ? x : size(a.strengths, dim)
    i1_filament_strs = replacefn(i1_filament_strs, 1)
    imax_filament_strs = replacefn(imax_filament_strs, 1)
    j1_filament_strs = replacefn(j1_filament_strs, 2)
    jmax_filament_strs = replacefn(jmax_filament_strs, 2)

    println(typeof(a))
    println(typeof(velocity_fn))
    pres = steady_pressures(a, velocity_fn, density,
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
        imax_filament_strs=imax_filament_strs,
        jmax_filament_strs=jmax_filament_strs,
        i1_filament_strs=i1_filament_strs,  
        j1_filament_strs=j1_filament_strs)
<<<<<<< HEAD
=======
    # We need forces back...
    mareas = areas(a.geometry)
    mforces = mareas .* pres
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
    mcents = centres(a.geometry)
    total_forces = sum(mforces)
    moment_contrib = map(
        x->cross(x[1]-measurement_centre, x[2]),
        zip(mcents, mforces))
<<<<<<< HEAD
    total_moment = sum(moment_contrib)
    # Calculate pressures... meh...
    mareas = areas(a.geometry)
    mnormals = normals(a.geometry)
    normal_forces = map(x->dot(x[1], x[2]), zip(mnormals, mforces))
    all_pressures = normal_forces ./ mareas
    return total_forces, total_moment, all_pressures
=======
    return total_forces, moment_contrib, pres
>>>>>>> 14419b807d5854ecfc5e86b34c55ef0698ce806e
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

    n = length(this.ring_vertices)
    @assert(length(state_vector) == n * 3, string(
        "State vector was incorrect length. Length was ", length(state_vect),
        " but should have been ", n * 3, "."))
    if any(isnan.(state_vect) || isinf.(state_vect))
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
    for i = 1 : length(this.ring_vertices)
        svtd[i*3 - 2 : i*3] = convert(Vector{Float64},
             induced_velocity(this, inducing_bodies))
    end
    return svtd
end

function vorticity_vector_length(this::VortexRingLattice)
    return length(this.strengths)
end

function vorticity_vector(this::VortexRingLattice)
    return vec(this.strengths)
end

function update_using_vorticity_vector!(
    this::VortexRingLattice, vort_vect::Vector{Float64})

    @assert(length(vort_vect) == vorticity_vector_length(this),
        string("Incorrect vorticity_vector length for VortexRingLattice.",
        " Expected a length of ", vorticity_vector_length(this), " but ",
        " was given a vector of length ", length(vort_vect)))
    for i = 1 : length(vort_vect)
        this.strengths[i] = vort_vect[i]
    end
    return
end

function vorticity_vector_velocity_influence(
    this :: VortexRingLattice,
    mes_pnt :: Vector3D)

    ret = zeros(3, vorticity_vector_length(this))
    indexes = [(i, j) for 
        i = 1 : size(this.strengths, 1), j = 1 : size(this.strengths, 2)]
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

#= CONVERSION FUNCTIONS ------------------------------------------------------=#
function Base.convert(::Type{Matrix{VortexRing}}, a::VortexRingLattice)
    ret = Matrix{VortexRing}(undef, size(a.strengths))
    g = a.geometry.coordinates
    s = a.strengths
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

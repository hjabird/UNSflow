#===============================================================================
    ThreeDVortexParticleSet

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDVortexParticleSet
    particles :: Vector{ThreeDVortexParticle}

    _reduction_factor_fn :: Function
    _vorticity_fraction_fn :: Function

    function ThreeDVortexParticleSet(
        particles :: Vector{ThreeDVortexParticle},
        reduction_factor_fn :: Function,
        vorticity_fraction_fn :: Function
        )

        _reduction_factor_fn = reduction_factor_fn
        _vorticity_fraction_fn = vorticity_fraction_fn
        return new(particles, reduction_factor_fn, vorticity_fraction_fn)
    end
end

function ThreeDVortexParticleSet(particles :: Vector{ThreeDVortexParticle})
    r, v = threed_winckelmans_kernels()
    return ThreeDVortexParticleSet(
        particles, r, v
    )
end

function ThreeDVortexParticleSet()
    r, v = threed_winckelmans_kernels()
    return ThreeDVortexParticleSet(
        Vector{ThreeDVortexParticle}(0), r, v
    )
end

function convert(
    ::Type{ThreeDVortexParticleSet},
    a::Vector{ThreeDVortexParticle})
    b = ThreeDVortexParticleSet(a)
    return b
end

function Base.length(a :: ThreeDVortexParticleSet)
    return length(a.particles)
end

function Base.size(a :: ThreeDVortexParticleSet)
    return Vector{Int64}([length(a)])
end

Base.eltype(::Type{ThreeDVortexParticleSet}) = ThreeDVortexParticle

Base.start(::ThreeDVortexParticleSet) = Int64(1)

function Base.next(a::ThreeDVortexParticleSet, state :: Int64)
    return (a.particles[state], state + 1)
end

function Base.done(a::ThreeDVortexParticleSet, state :: Int64)
    return state > length(a)
end

function Base.endof(a::ThreeDVortexParticleSet)
    return length(a.particles)
end

function Base.getindex(a::ThreeDVortexParticleSet, i :: Int)
    return a.particles[i]
end

function Base.getindex(a::ThreeDVortexParticleSet, I)
    return [a.particles[s] for s in I]
end

function Base.setindex!(a::ThreeDVortexParticleSet,
    val::ThreeDVortexParticle, i :: Int)
    a[i] = val
end

function Base.:+(a::ThreeDVortexParticleSet, b::ThreeDVortexParticleSet)
    if(a._reduction_factor_fn != b._reduction_factor_fn)
        error(string("Both vortex particle sets must be using the same ",
            "reduction functions to be able to merge."))
    end
    if(a._vorticity_fraction_fn != b._vorticity_fraction_fn)
        error(string("Both vortex particle sets must be using the same ",
            "vorticity fraction functions to be able to merge."))
    end
    c = a
    c.particles = vcat(a.particles, b.particles)
    return c
end

function Base.:+(a::ThreeDVortexParticleSet, b::Vector{ThreeDVortexParticle})
    c = a
    c.particles = vcat(a.particles, b)
    return c
end

function Base.:+(a::Vector{ThreeDVortexParticle}, b::ThreeDVortexParticleSet)
    c = b
    c.particles = vcat(a, b.particles)
    return c
end
#= END ThreeDVortexParticleSet -----------------------------------------------=#

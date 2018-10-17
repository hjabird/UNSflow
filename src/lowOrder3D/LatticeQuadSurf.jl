#===============================================================================
    LatticeQuadSurf.jl

    A surface structured lattice of quad rings.

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

mutable struct LatticeQuadSurf <: DiscreteGeometry3D
    coordinates :: Matrix{Vector3D}

    function LatticeQuadSurf(a::Matrix{Vector3D})
        new(a)
    end
    function LatticeQuadSurf(isize::Int, jsize::Int)
        c = Matrix{Vector3D}(undef, (isize, jsize))
        new(c)     
    end
end

# Return map a local coordinate to a point in space
function evaluate(a::LatticeQuadSurf, local_coord::Vector{T},
    surfaceinterpolating::Bool=false) where T <: Real
    @assert(length(position) == 2, "LatticeQuadSurf is two dimensionsal.")
    if !surfaceinterpolating
        @assert(all(abs.(local_coord - round.(local_coord) .<= eps(T)),
            string("The LatticeQuadSurf is only defined on the lattice ",
                "section. IE, it has holes in it compared to BilinearQuadSurf",
                ". Either use BilinearQuadSurf, or - if you must - call the",
                " evaluate function as evaluate(::LatticeQuadSurf, local_coord",
                "::Vector{T}, true). Tried to interpolate at ", local_coord)))
    end
    xf = Int64(local_coord[1])
    yf = Int64(local_coord[2])
    v = a.coordinates
    g = [v[xf+1, yf], v[xf+1, yf+1], v[xf, yf+1], v[xf, yf]]
    x = (local_coord[1] % 1) * 2 - 1
    y = (local_coord[2] % 1) * 2 - 1
    t1 = g[1] * (x-1)(y-1) / 4
    t2 = - g[2] * (x+1)(y-1) / 4
    t3 = g[3] * (x+1)(y+1) / 4
    t4 = - g[4] * (x-1)(y+1) / 4
    return t1 + t2 + t3 + t4
end

# Return a vector of coordinates the interpolation points of the geometry
function coords(a::LatticeQuadSurf)
    return vec(a.coordinates)
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::LatticeQuadSurf)
    return 2
end

# Get the number of control points of an object
function number_of_control_points(a::LatticeQuadSurf)
    return length(a.coordinates)
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::LatticeQuadSurf, position::Vector{T},
    surfaceinterpolating::Bool=false) where T <: Real
    @assert(length(position) == 2, "LatticeQuadSurf is two dimensionsal.")
    if !surfaceinterpolating
        return all(abs.(local_coord - round.(local_coord) .<= eps(T)))
    end
    return (1 <= position[1] <= size(a.coordinates, 1)) &&
        (1 <= position[2] <=  size(a.coordinates, 2))
end

function Base.size(a::LatticeQuadSurf)
    i, j = size(a.coordinates)
    return (i-1, j-1)
end

function Base.convert(::Type{PolyLine2}, a::LatticeQuadSurf,
    xidx::Int, yidx::Int)

    @assert(xidx < size(a.coordinates, 1), string("xidx of requested PolyLine2",
        " was too high. Requested x index of ", xidx, 
        " which should have been less than ", size(a.coordinates, 1), "."))
    @assert(yidx < size(a.coordinates, 2), string("yidx of requested PolyLine2",
        " was too high. Requested y index of ", yidx, 
        " which should have been less than ", size(a.coordinates, 2), "."))
    v = a.coordinates
    g = [v[xidx+1, yidx], v[xidx+1, yidx+1], v[xidx, yidx+1], v[xidx, yidx],
         v[xidx+1, yidx]]
    return PolyLine2(g)
end

function Base.convert(::Type{Matrix{PolyLine2}}, a::LatticeQuadSurf)
    indexes = [(i, j) for i = 1 : size(a.coordinates, 1) - 1,
                          j = 1 : size(a.coordinates, 2) - 1]
    ret = map(x->convert(PolyLine2, a, x[1], x[2]), indexes)
    return ret
end

function Base.convert(::Type{Vector{PolyLine2}}, a::LatticeQuadSurf)
    return vec(convert(Matrix{PolyLine2}, a))
end

function Base.convert(::Type{Vector{Line2}}, a::LatticeQuadSurf,
    xidx::Int, yidx::Int)
    return convert(Vector{Line2}, convert(PolyLine2, a, xidx, yidx))
end

function Base.convert(::Type{Vector{Line2}}, a::LatticeQuadSurf)
    ci, cj = size(a.coordinates)
    ci -= 1
    cj -= 1
    n = ci * cj * 2 + ci + cj
    ret = Vector{Line2}(undef, n)
    # The top and left sides of each "box"
    for i = 1 : ci 
        for j = 1 : cj
            lind = ((i - 1) * cj + j) * 2 - 1
            line2s = convert(Vector{Line2}, a, i, j)
            ret[lind] = line2s[4]   # Vertical, left
            ret[lind + 1] = line2s[3]   # Horizontal, top
        end
    end
    offset = ci * cj * 2
    # The right edge
    for i = 1 : ci
        line2s = convert(Vector{Line2}, a, i, cj)
        ret[offset + i] = line2s[2]
    end
    offset += ci
    # The bottom edge
    for j = 1 : cj
        line2s = convert(Vector{Line2}, a, ci, j)
        ret[offset + j] = line2s[1]
    end
    return ret
end

function Base.convert(::Type{BilinearQuadSurf}, a::LatticeQuadSurf)
    return BilinearQuadSurf(a.coordinates)
end

# DEVELOPERS:
#   This function is out of place - it ought to live in BilinearQuadSurf.jl,
#   but due to the unique way Julia development ignores issues (see issue 
#   269, open since 2011) this can't live where it ought to. 
function Base.convert(::Type{LatticeQuadSurf}, a::BilinearQuadSurf)
    return LatticeQuadSurf(a.coordinates)
end

function normals(a::LatticeQuadSurf)
    v = a.coordinates
    ret = Matrix{Vector3D}(undef, size(a))
    for i = 1 : size(ret, 1)
        for j = 1 : size(ret, 2)
            g = [v[i+1, j], v[i+1, j+1], v[i, j+1], v[i, j]]
            t1 = 0.25 * (g[4] + g[3] - g[2] - g[1])
            t2 = 0.25 * (g[2] + g[3] - g[1] - g[4])
            ret[i, j] = unit(cross(t1, t2))
        end
    end
    return ret
end

function centres(a::LatticeQuadSurf)   
    v = a.coordinates
    ret = Matrix{Vector3D}(undef, size(a))
    for i = 1 : size(ret, 1)
        for j = 1 : size(ret, 2)
            g = [v[i+1, j], v[i+1, j+1], v[i, j+1], v[i, j]]
            ret[i, j] = sum(g) / 4
        end
    end
    return ret
end

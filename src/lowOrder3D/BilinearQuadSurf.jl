#===============================================================================
    BilinearQuadSurf.jl

    A surface made up of a structured grid of BilinearQuads

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

mutable struct BilinearQuadSurf <: DiscreteGeometry3D
    coordinates :: Matrix{Vector3D}

    function BilinearQuadSurf(a::Matrix{Vector3D})
        new(a)
    end
    function BilinearQuadSurf(isize::Int, jsize::Int)
        c = Matrix{Vector3D}(undef, (isize, jsize))
        new(c)     
    end
end

# Return map a local coordinate to a point in space
function evaluate(a::BilinearQuadSurf, local_coord::Vector{T}) where T <: Real
    @assert(length(position) == 2, "BilinearQuadSurf is two dimensionsal.")
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
function coords(a::BilinearQuadSurf)
    return vec(a.coordinates)
end

# Get the number of dimensions that the space operates in (ie, line->1, surf->2)
function local_dimensions(a::BilinearQuadSurf)
    return 2
end

# Get the number of control points of an object
function number_of_control_points(a::BilinearQuadSurf)
    return length(a.coordinates)
end

# Test if a point is in the bounds defined by the object
function in_bounds(a::BilinearQuadSurf, position::Vector{T}) where T <: Real
    @assert(length(position) == 2, "BilinearQuadSurf is two dimensionsal.")
    return (1 <= position[1] <= size(a.coordinates, 1)) &&
        (1 <= position[2] <=  size(a.coordinates, 2))
end

function Base.size(a::BilinearQuadSurf)
    i, j = size(a.coordinates)
    return (i-1, j-1)
end

function Base.convert(::Type{BilinearQuad}, a::BilinearQuadSurf,
    xidx::Int, yidx::Int)

    @assert(xidx < size(a.coordinates, 1), string("xidx of requested Bilinear",
        "Quad was too high. Requested x index of ", xidx, 
        " which should have been less than ", size(a.coordinates, 1), "."))
    @assert(yidx < size(a.coordinates, 2), string("yidx of requested Bilinear",
        "Quad was too high. Requested y index of ", yidx, 
        " which should have been less than ", size(a.coordinates, 2), "."))
    v = a.coordinates
    g = [v[xidx+1, yidx], v[xidx+1, yidx+1], v[xidx, yidx+1], v[xidx, yidx]]
    return BilinearQuad(g)
end

function Base.convert(::Type{Matrix{BilinearQuad}}, a::BilinearQuadSurf)
    nx, ny = size(a.coordinates)
    ret = map(x->convert(BilinearQuad, a, x[1], x[2]),
        [(i, j) for i = 1 : size(a.coordinates, 1) - 1, 
            j = 1 : size(a.coordinates, 2) - 1])
    return ret
end

function Base.convert(::Type{Vector{BilinearQuad}}, a::BilinearQuadSurf)
    return vec(convert(Matrix{BilinearQuad}, a))
end

function normals(a::BilinearQuadSurf)
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

function centres(a::BilinearQuadSurf)   
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
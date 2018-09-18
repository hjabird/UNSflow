#===============================================================================
    ThreeDVortexRing

    Representation of a singular vortex filament ring.

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
type ThreeDVortexRing
    c1 :: ThreeDVector
    c2 :: ThreeDVector
    c3 :: ThreeDVector
    c4 :: ThreeDVector
    strength :: Float64

    function ThreeDVortexRing(
        corner1 :: ThreeDVector,
        corner2 :: ThreeDVector,
        corner3 :: ThreeDVector,
        corner4 :: ThreeDVector,
        strength :: Float64
        )
        return new(corner1, corner2, corner3, corne4, strength)
    end
end

function ThreeDVortexRing(
    corner1 :: ThreeDVector,
    corner2 :: ThreeDVector,
    corner3 :: ThreeDVector,
    corner4 :: ThreeDVector
    )
    return ThreeDVortexRing(corner1, corner2, corner3, corner4, 0.0)
end

function ThreeDVortexRing()
    c1 = ThreeDVector()
    c2 = ThreeDVector()
    c3 = ThreeDVector()
    c4 = ThreeDVector()
    return ThreeDVortexRing(c1, c2, c3, c4, 0.0)
end

function convert(::Vector{ThreeDStraightVortexFilament}, a::ThreeDVortexRing)
    b = Vector{ThreeDStraightVortexFilament}(4)
    b[1] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    b[2] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    b[3] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    b[4] = ThreeDStraightVortexFilament(a.c1, a.c2, a.strength)
    return b
end
#= END ThreeDVortexRing ------------------------------------------------------=#

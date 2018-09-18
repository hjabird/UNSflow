#===============================================================================
    SplineChordShape.jl

    A type to define a wing chord line geometry in 2D. Defined in [-1, 1] for
    leading edge to trailing edge. The y location is 0 at both LE at TE.

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

type SplineChordShape
    camber_line :: Dierckx.Spline1D # Y locations for x in [-1, 1]

    function SplineChordShape(
        camber_line :: Spline1D
        )
        @assert(abs(camber_line(-1)) <= eps{Float64}(),
            "The leading edge (-1) must be at y=0")
        @assert(abs(camber_line(1)) <= eps{Float64}(),
            "The trailing edge (1) must be at y=0")
        new(camber_line)
    end
end

# Default constructor returns a flat chord geometry.
function SplineChordShape()
    xs = [-1., 1.]
    ys = [0., 0.]
    camber_line = Spline1D(xs, ys, k=1)
    return ChordShape(camber_line)
end

function chord_length(chord::SplineChordShape)
    return 2.
end

function y_location(chord::SplineChordShape, x_position::Number)
    @assert(-1 <= x_position <= 1, "The chord is defined in [-1, 1]")
    return Dierckx.eval(chord.camber_line, x_position)
end

function dydx(chord::SplineChordShape, x_position::Number)
    @assert(-1 <= x_position <= 1, "The chord is defined in [-1, 1]")
    return Dierckx.derivative(chord.camber_line, x_position)
end

#= END WingChordSection ------------------------------------------------------=#

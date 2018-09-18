#===============================================================================
    ThreeDVorticityBody.jl

    Represents some kind of vorticity within the flow domain. Does not define
    goemetry or any other such thing.

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

abstract type ThreeDVorticityBody
end

# Expected interface for 3DVBody:

# Constructor is defined by the concrete type.

# centre(this::3DVBody)
#   returns something somehow representative of the bodies' centre location
#   as a ThreeDVector. Whether this is geometrically central, or vorticity
#   weighted does not matter so long as it is considered by the programmer along
#   with the effective_radius definition.

# effective_radius(this::3DBVody)
#   returns a Real representative of the effective radius from the centre
#   of a sphere that contains the 3D vorticity body's vorticy.

# vorticity(this::3DVBody)
#   returns the integral of vorticity within the body as ThreeDVector.

# induced_velocity(this::3DVBody, measurement_point::ThreeDVector)
#   returns the velocity induced by the body at measurement_point as a
#   ThreeDVector.

# induced_velocity_curl(this::3DVBody, measurement_point::ThreeDVector)
#   returns the curl in the velocity induced by the body at measurement_point as
#   a 3 by 3 matrix, by which vorticity at that point can be multiplied.

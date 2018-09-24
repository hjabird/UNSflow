#===============================================================================
    RedistributionScheme.jl

    Redistribution schemes for vortex particles as described in Winckelmans,
    C. R. Physique 6 (2005).

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

function zeroeth_order_redistribution_scheme(U::Real)
    if U < 0.5
        return 1.
    else
        return 0.
    end
end

function first_order_redistribution_scheme(U::Real)
    if 0 <= U <= 1
        return 1 - U
    else
        return 0.
    end
end

function second_order_redistribution_scheme(U::Real)
    if U < 0 || U >= 1.5
        return 0.
    elseif U < 0.5
        return 1 - U^2
    else
        return 0.5 * (1 - U) * (2 - U)
    end
end

function third_order_redistribution_scheme(U::Real)
    if U < 0 || U > 2
        return 0.
    elseif U <= 1
        return 0.5 * (1 - U^2) * (2 - U)
    else
        return (1/6) * (1 - U) * (2 - U) * (3 - U)
    end
end

function m4_redistribution_scheme(U::Real)
    if U < 0 || U > 2
        return 0.
    elseif U <= 1
        return 1 - (5/2)*U^2 + (3/2)*U^3
    else
        return 0.5 * (1 - U) * (2 - U)^2
    end
end

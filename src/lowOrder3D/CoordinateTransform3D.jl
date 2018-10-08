#===============================================================================
    CoordinateTransform3D.jl

    A time variant transplation of a points in space.

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

mutable struct CoordinateTransform3D
    func :: Function
    time :: Float64
    derivative_method :: Function

    function CoordinateTransform3D(
        func :: Function,
        time :: Float64,
        derivative_method :: Function)
        return new(func, time, derivative)
    end
end

function CoordinateTransform3D(
    func :: Function,
    time :: Float64)
    t = CoordinateTransform3D(
        func :: Function,
        time :: Float64,
        x->x)
    set_central_difference!(t, 1e-6)
    return t
end

function CoordinateTransform3D(
    func :: Function)
    t = CoordinateTransform3D(
        func :: Function,
        0.0)
    return t
end

function CoordinateTransform3D()
    t = CoordinateTransform3D((x,t)->x)
    return t
end

function (t::CoordinateTransform3D)(x::Vector3D)
    return t.func(x, t.time)
end

function evaluate(t::CoordinateTransform3D, x::Vector3D)
    return t(x)
end

function derivative(t::CoordinateTransform3D, x::Vector3D)
    return t.derivative_method(x, t.time)
end

function set_central_difference!(ct::CoordinateTransform3D, dt::Float64)
    ct.derivative_method =  (x::Vector3D, t::Float64)->
        (ct.func(x, t + dt / 2) - ct.func(x, t - dt / 2)) / dt
end

function set_backward_difference!(ct::CoordinateTransform3D, dt::Float64)
    ct.derivative_method = (x::Vector3D, t::Float64)->
                                (ct.func(x, t) - ct.func(x, t - dt)) / dt
end

function set_derivative_expression!(
    t::CoordinateTransform3D,
    expr :: Function)
    t.derivative_method = expr
end

function increment!(t::CoordinateTransform3D, dt::Float64)
    return t.time += dt
end

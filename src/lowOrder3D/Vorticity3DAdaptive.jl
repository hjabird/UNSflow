#===============================================================================
    Vorticity3DAdaptive.jl

    Vortex bodies that are in some way adaptive. Primarily this is for
    implementing vortex particle redistribution, but is also applicable to
    vortex filament spline approximations etc.

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
include("Vorticity3DCollector.jl")

abstract type Vorticity3DAdaptive <: Vorticity3DCollector
end

# This should be reimplemented by each adaptive container to reflect its own
# adaptive nature. If -oddly- you wanted to knowthing, this should be
# reimplemented to do nothing.
function adaptive_update!(a::Vorticity3DAdaptive)
    error(string("The function adaptive_update(a::", typeof(a),
        ") has not been reimpleneted for that method as it should have!"))
end

function Base.push!(a::Vorticity3DAdaptive, to_be_added::Vorticity3D)
    error("Adaptive collectors require additional information to add
        elements. Either look in the source for the method or, in the Julia
        REPL, try ?methodswith(<type_to_examine>).")
end

function Base.append!(
    a::Vorticity3DAdaptive,
    iterable_of_things_to_be_appended)
    error("Adaptive collectors require additional information to add
        elements. Either look in the source for the method or, in the Julia
        REPL, try ?methodswith(<type_to_examine>).")
end


function Base.pop!(a::Vorticity3DAdaptive)
    error("Popping a Vorticity3DAdaptive is not allowed by default.
    Add a concrete instance specific method if one is not already provided.")
end

function Base.isempty(a::Vorticity3DAdaptive)
    return isempty(a.children)
end

function Base.empty!(a::Vorticity3DAdaptive)
    error("Emptying a Vorticity3DAdaptive is not allowed. If you want rid
    of its contents, just delete it instead. Alternately, if it makes sense,
    add a concrete instance specific method if one is not already provided.")
end

#= END Vorticity3DAdaptive ---------------------------------------------------=#

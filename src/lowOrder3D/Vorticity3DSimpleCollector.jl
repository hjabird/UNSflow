#===============================================================================
    Vorticity3DSimpleCollector.jl

    A container for Vorticity3D(ies). As simple a container as is
    possible.

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

# All methods are inherited from the Vorticity3DCollector supertype.

mutable struct Vorticity3DSimpleCollector <: Vorticity3DCollector
    children :: Vector{Vorticity3D}

    function Vorticity3DSimpleCollector()
        new(Vector{Vorticity3D}())
    end
    function Vorticity3DSimpleCollector(child::Vorticity3DCollector)
        new(child)
    end
    function Vorticity3DSimpleCollector(children...)
        a = Vorticity3DSimpleCollector()
        for child in children
            if typeof(child) <: Vorticity3D
                push!(a.children, child)
            else
                try
                    for childchild in child
                        push!(a.children, childchild)
                    end
                catch
                    error("Children must be Vorticity3D objects or
                        iterable containers of Vorticity3D objects.")
                end
            end
        end
        return a
    end
    function Vorticity3DSimpleCollector(children_iterable)
        c = Vector{Vorticity3D}()
        for child in children_iterable
            push!(c, child)
        end
        new(c)
    end
end

#= END Vorticity3DSimpleCollector --------------------------------------------=#

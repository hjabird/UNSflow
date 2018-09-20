#===============================================================================
    ThreeDVorticitySimpleCollector.jl

    A container for ThreeDVorticityBody(ies). As simple a container as is
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
include("ThreeDVorticityCollector.jl")

# All methods are inherited from the ThreeDVorticityCollector supertype.

type ThreeDVorticitySimpleCollector <: ThreeDVorticityCollector
    children :: Vector{ThreeDVorticity}

    function ThreeDVorticitySimpleCollector()
        new(Vector{ThreeDVorticity}())
    end
    function ThreeDVorticitySimpleCollector(child::ThreeDVorticityCollector)
        new(child)
    end
    function ThreeDVorticitySimpleCollector(children...)
        a = ThreeDVorticitySimpleCollector()
        for child in children
            if typeof(child) <: ThreeDVorticity
                push!(a.children, child)
            else
                try
                    for childchild in child
                        push!(a.children, childchild)
                    end
                catch
                    error("Children must be ThreeDVorticity objects or
                        iterable containers of ThreeDVorticity objects.")
                end
            end
        end
        return a
    end
    function ThreeDVorticitySimpleCollector(children_iterable)
        c = Vector{ThreeDVorticity}()
        for child in children_iterable
            push!(c, child)
        end
        new(c)
    end
end

#= END ThreeDVorticitySimpleCollector -------------------------------------=#

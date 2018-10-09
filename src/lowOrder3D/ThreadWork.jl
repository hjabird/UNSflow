#===============================================================================
    ThreadWork.jl

    An object to help a multi-threaded setup navigate the trees generated
    in Vorticity3D to acceleration computation.

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

mutable struct ThreadWorkload
    maxidx :: Int64     # This should stay constant
    nextindex :: Threads.Atomic{Int64}  # Needs lock to avoid read before write etc.
    children :: Vector{ThreadWorkload}
    child_lock :: Threads.SpinLock  # Lock readwrite access to children vect
    worker_count :: Threads.Atomic{Int64}
    object :: Any   # The object in the tree that the workload is paired with.

    # You MUST enter a workload upon its creation!
    function ThreadWorkload(
        object :: Any
        )

        new(length(object), Threads.Atomic{Int64}(1), Vector{ThreadWorkload}(),
            Threads.SpinLock(), Threads.Atomic{Int64}(0), object)
    end
end

mutable struct ThreadState
    workloads :: Vector{ThreadWorkload}
end

function next_idx!(a::ThreadWorkload)
    w = Threads.atomic_add!(a.nextindex, 1)
    if w > a.maxidx
        w = -1
    end
    return w
end

function add_child!(a::ThreadState, child_workload::ThreadWorkload)
    Threads.lock(a.workloads[end].child_lock)
    push!(a.workloads[end].children, child_workload)
    Threads.unlock(a.workloads[end].child_lock)
    enter_TheadWorkload!(a, child_workload)
    return
end

function enter_TheadWorkload!(a::ThreadState, b::ThreadWorkload)
    push!(a.workloads, b)
    Threads.atomic_add!(b.worker_count, 1)
end

function exit_ThreadWorkload!(a::ThreadState)
    current = a.workloads[end]
    pop!(a.workloads)
    Threads.atomic_sub!(current.worker_count, 1)
    if current.worker_count.value == 0
        Threads.lock(a.workloads[end].child_lock)
        to_remove = findall(x->current===x, a.workloads[end].children)
        deleteat!(a.workloads[end].children, to_remove[1])
        Threads.unlock(a.workloads[end].child_lock)
        @assert(length(to_remove) == 1)
    end
end

function apply_to_tree(func::Function, root::Any, iterable_supertype)
    local worker_func
    function worker_func()
        local tstate = ts[Threads.threadid()]
        while length(tstate.workloads) > 0
            workld = tstate.workloads[end]
            idx = next_idx!(workld)
            if idx < 0
                if length(workld.children) > 0
                    break
                    enter_TheadWorkload!(tstate, workld.children[1])
                else
                    exit_ThreadWorkload!(tstate)
                end
            elseif typeof(workld.object[idx]) <: iterable_supertype
                add_child!(tstate, ThreadWorkload(workld.object[idx]))
            else
                func(workld.object[idx])
            end
        end
    end
    ts = Vector{ThreadState}(undef, Threads.nthreads())
    for i = 1 : Threads.nthreads()
        startstate = ThreadState([ThreadWorkload(root)])
        ts[i] = startstate
    end
    ccall(:jl_threading_run, Ref{Cvoid}, (Any,), worker_func)
    return
end

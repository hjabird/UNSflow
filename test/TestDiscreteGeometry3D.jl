#===============================================================================
    TestDiscreteGeometry3D.jl

    Tests for TestDiscreteGeometry3D stuff

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

push!(LOAD_PATH, "../src/")
import UNSflow
import Test

let
    bq1 = UNSflow.BilinearQuad(UNSflow.Vector3D(0,0,0), 
        UNSflow.Vector3D(1,0,0),  UNSflow.Vector3D(1,1,0), 
        UNSflow.Vector3D(0,1,0))
    bq1c = deepcopy(bq1)
    bq2 =  UNSflow.BilinearQuad( UNSflow.Vector3D(0,0,0),  
        UNSflow.Vector3D(-1,0,0), UNSflow.Vector3D(-1,-1,0),  
        UNSflow.Vector3D(0,-1,0))

    Test.@test bq1 == bq1
    Test.@test bq1 == bq1c
    Test.@test bq1 != bq2
    Test.@test isequal(bq1, bq1c)
    Test.@test isequal(bq1, bq2) == false
    Test.@test hash(bq1) == hash(bq1)
    Test.@test hash(bq1) == hash(bq1c)
    Test.@test hash(bq1) != hash(bq2)
end


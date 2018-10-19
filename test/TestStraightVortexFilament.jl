#===============================================================================
    TestStraightVortexFilament.jl

    Tests for StraightVortexFilament

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
    extv0 = UNSflow.Vector3D(0, 0, 0)
    extvx1 = UNSflow.Vector3D(1, 0, 0)
    extvy1 = UNSflow.Vector3D(0, 1, 0)
    extvz1 = UNSflow.Vector3D(0, 0, 1)


    # Vorticity calculation
    let
        vf1 = UNSflow.StraightVortexFilament(extv0, extvx1, 1)
        Test.@test UNSflow.vorticity(vf1) == extvx1
        vf2 = UNSflow.StraightVortexFilament(extv0, extvx1, 2)
        Test.@test UNSflow.vorticity(vf2) == extvx1 * 2
        vf3 = UNSflow.StraightVortexFilament(extv0, extvy1, 1)
        Test.@test UNSflow.vorticity(vf3) == extvy1
    end

    # Induced velocity
    let
        vf1 = UNSflow.StraightVortexFilament(extv0, extvx1, 1)
        Test.@test UNSflow.induced_velocity(vf1, extv0) == extv0
        Test.@test UNSflow.induced_velocity(vf1, extvx1) == extv0
        Test.@test UNSflow.induced_velocity(vf1, extvx1 * 0.5) == extv0
        Test.@test UNSflow.induced_velocity(vf1, extvy1 * 1e-15) == extv0
        Test.@test UNSflow.induced_velocity(vf1, extvy1 * 0.5) != extv0
    end

    # Geometry
    let
        vf1 = UNSflow.StraightVortexFilament(extv0, extvx1, 1)
        Test.@test UNSflow.centre(vf1) == extvx1 / 2
        Test.@test UNSflow.effective_radius(vf1) == 0.5
    end

    # steady force function.
    let
        vel_fn_ux = x->extvx1
        vel_fn_uz = x->extvz1
        Test.@test UNSflow.steady_force(UNSflow.StraightVortexFilament,
            extv0, extvx1, 1., vel_fn_ux) == extv0
        Test.@test UNSflow.steady_force(UNSflow.StraightVortexFilament,
            extv0, extvx1, 1., vel_fn_uz) == extvy1
        Test.@test UNSflow.steady_force(UNSflow.StraightVortexFilament,
            extv0, extvx1, 2., vel_fn_uz) == extvy1 * 2
        Test.@test UNSflow.steady_force(UNSflow.StraightVortexFilament,
            extv0, extvx1, 1., vel_fn_uz, 2.) == extvy1 * 2
        Test.@test UNSflow.steady_force(UNSflow.StraightVortexFilament,
            extv0, extvx1, 1., vel_fn_uz, 1, 2) == extvy1
    end
end # let

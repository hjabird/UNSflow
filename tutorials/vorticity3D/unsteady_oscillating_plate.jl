#===============================================================================
    unsteady_vortex_lattice_method.jl

    An example of using the vortex lattice method.

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
push!(LOAD_PATH, "../../src/")
import UNSflow

let
# Define a rectangular plate:
surf_fn = x->UNSflow.Vector3D(x[1], x[2], 0)
disc = map(x->cbrt(x), -1:0.2:1)
plate_surf = UNSflow.EquationSurf(surf_fn)
plate_geom = UNSflow.discretise(plate_surf, UNSflow.BilinearQuadSurf, 
                                                        disc, disc)
# Define some kinematics:
# kinem = UNSflow.CoordinateTransform3D((x,t)->
#    UNSflow.Vector3D(0, 0, 0.25*cos(t*pi)) + x)
kinem = UNSflow.CoordinateTransform3D((x,t)->
    UNSflow.Vector3D(0, 0, 0.1 * t * t) + x)
# We use deepcopy since we need a "baseline" to apply the transform to.
plate = UNSflow.VortexRingLattice(deepcopy(plate_geom))
dt = 0.05

# When we solve the problem, its a Neumann problem.
problem = UNSflow.NeumannProblem3D()
push!(problem.variable_vorticity, plate)
max=200
# We want the old vorticity for when we work out the unsteady loads.
old_vort_vect = UNSflow.vorticity_vector(plate)
for i = 1 : max
    # We have to redo the boundary conditions on each step since they change.
    UNSflow.clear_bcs!(problem)
    wc_coords = plate.geometry
    wc_coords = map(x->kinem(x), UNSflow.coords(plate_geom))
    wcentres = UNSflow.centres(plate.geometry)
    wnormals = UNSflow.normals(plate.geometry)
    wvel = map(x->UNSflow.derivative(kinem, x), UNSflow.centres(plate_geom))
    wnvel = map(x->UNSflow.dot(UNSflow.unit(x[1]), x[2]), zip(wnormals, wvel))
    UNSflow.add_bc!(problem, wcentres, wnormals, wnvel)

    # Solve the problem
    UNSflow.solve!(problem)

    # Total loads = (psuedo)steady load + unsteady loads
    println("Step ", i)
    ind_vel_fn = x -> UNSflow.induced_velocity(plate, x) - 
        UNSflow.derivative(kinem, x)
    sforce, smoment, spressure = UNSflow.steady_loads(
        problem.variable_vorticity[1], ind_vel_fn)
    usforce, usmoment, uspressure = UNSflow.unsteady_loads(
        plate, old_vort_vect, dt)
    new_vort_vect = UNSflow.vorticity_vector(plate)
    old_vort_vect = new_vort_vect
        
    println("Forces steady: ", sforce)
    println("Forces total: ", usforce + sforce)
    
    # We'll use the MeshDataLinker to add the steady and and unsteady
    # pressures to the problem.
    extra_data = UNSflow.MeshDataLinker()
    UNSflow.add_celldata!(extra_data, problem.variable_vorticity[1], 
        "SteadyPressure", spressure)
    UNSflow.add_celldata!(extra_data, problem.variable_vorticity[1], 
        "UnsteadyPressure", uspressure)
    UNSflow.add_celldata!(extra_data, problem.variable_vorticity[1], 
        "TotalPressure", uspressure + spressure)
    UNSflow.add_celldata!(extra_data, problem.variable_vorticity[1], 
        "NormalVelBC", map(x->UNSflow.unit(x[1])*x[2], zip(wnormals, wnvel)))
    mesh = UNSflow.UnstructuredMesh()
    push!(mesh, UNSflow.vorticities(problem))
    UNSflow.add_data!(mesh, extra_data)
    UNSflow.to_vtk_file(mesh, string("output/Unsteady_oscillating_plate_", i))

    # And its so easy to forget this kind of thing...
    UNSflow.increment!(kinem, dt)
    plate.geometry.coordinates = map(x->kinem(x), plate_geom.coordinates)
end

end #let

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

println("GOING...")
let
# STEP 1 -----------------------------------------------------------------------
aoa = 4 * pi / 180
aspect_ratio = 3
# Rectangular:
surf_fn = x->UNSflow.Vector3D(x[1]/2, x[2]*aspect_ratio/2, 0)
# Or perhaps elliptic:
#surf_fn = x->UNSflow.Vector3D(sqrt(1- x[2]^2) * x[1]/2, x[2]*aspect_ratio/2, 0)
# Apply the angle of attack:
surf_fn_aoa = x->UNSflow.rotate_about_y(surf_fn(x), aoa)
# The generate the surface and discretise the surface:
wing_surf = UNSflow.EquationSurf(surf_fn_aoa)
spanwise_disc = map(x->cbrt(x), -1:0.2:1)
chordwise_disc = map(x->x, -1 : 0.125: 1)
wing_geom = UNSflow.discretise(wing_surf, UNSflow.BilinearQuadSurf, 
                                                chordwise_disc, spanwise_disc)
# Free stream / kinematics:
free_stream = UNSflow.FreeStream3D(0, 0, 0)
kinem = UNSflow.CoordinateTransform3D((x,t)->UNSflow.Vector3D(-t, 0, 0.0*cos(t*pi/2)) + x)
kutta_wake = UNSflow.VortexRingLattice(1, length(spanwise_disc)-1)
wing = UNSflow.VortexRingLattice(deepcopy(wing_geom))
dt = 0.05

wing_wake_kutta_idxs = UNSflow.extract_vorticity_vector_indexes(
    wing, [size(wing.geometry, 1)], collect(1 :size(wing.geometry, 2)))
wing_wake_kutta_mtrx = zeros(   length(kutta_wake.vorticity), 
                                length(wing.vorticity))
for idx = zip(collect(1:20), wing_wake_kutta_idxs)
    wing_wake_kutta_mtrx[idx[1], idx[2]] = 1
end
linked_wake = UNSflow.LinkedVorticity3DChild(kutta_wake)
linked_wing = UNSflow.LinkedVorticity3DParent(
    wing, linked_wake, wing_wake_kutta_mtrx)
linked_wake.geometry.coordinates[1,:] = 
    map(x->kinem(x) + UNSflow.induced_velocity(free_stream, x) * dt, 
    wing.geometry.coordinates[end,:])
linked_wake.geometry.coordinates[2,:] = 
    map(x->kinem(x; time=-dt) , wing.geometry.coordinates[end,:])
free_wake = UNSflow.VortexRingLattice(0, length(spanwise_disc) - 1)
free_wake.geometry.coordinates[1, :] = linked_wake.geometry.coordinates[end, :]

problem = UNSflow.NeumannProblem3D()
push!(problem.variable_vorticity, linked_wing)
push!(problem.variable_vorticity, linked_wake)
push!(problem.invariant_vorticity, free_wake)
push!(problem.invariant_vorticity, free_stream)
max=150
old_vort_vect = UNSflow.vorticity_vector(problem.variable_vorticity[1])
for i = 1 : max
    # Set up Neumann bc on wing surface:
    UNSflow.clear_bcs!(problem)
    wc_coords = wing.geometry
    wc_coords = map(x->kinem(x), UNSflow.coords(wing_geom))
    wcentres = UNSflow.centres(wing.geometry)
    wnormals = UNSflow.normals(wing.geometry)
    wvel = map(x->UNSflow.derivative(kinem, x), UNSflow.centres(wing_geom))
    wnvel = map(x->UNSflow.dot(UNSflow.unit(x[1]), x[2]), zip(wnormals, wvel))
    UNSflow.add_bc!(problem, wcentres, wnormals, wnvel)

    # Solve and convect
    UNSflow.solve!(problem)

    ind_vel_fn = x -> UNSflow.induced_velocity(UNSflow.vorticities(problem), x) - 
        UNSflow.derivative(kinem, x)
    sforce, smoment, spressure = UNSflow.steady_loads(
        UNSflow.source_object(problem.variable_vorticity[1]), ind_vel_fn; 
        # The wing vortex lattice is overlapped by that of the wake, so we need to 
        # correct for that - imax indicates the maximum i index of the lattice.
        imax_filament_strs=-problem.variable_vorticity[2].vorticity[1,:])
    usforce, usmoment, uspressure = UNSflow.unsteady_loads(
        UNSflow.source_object(problem.variable_vorticity[1]), old_vort_vect, dt)
    new_vort_vect = UNSflow.vorticity_vector(problem.variable_vorticity[1])
    old_vort_vect = new_vort_vect
    wake_vels = map(x->UNSflow.induced_velocity(UNSflow.vorticities(problem), x), UNSflow.coords(free_wake.geometry))
    
    println("Forces steady: ", sforce)
    println("Forces total: ", usforce + sforce)

    println("Step ", i)
    
    extra_data = UNSflow.MeshDataLinker()
    UNSflow.add_celldata!(extra_data, problem.variable_vorticity[1], 
        "SteadyPressure", spressure)
    UNSflow.add_celldata!(extra_data, problem.variable_vorticity[1], 
        "TotalPressure", spressure + uspressure)
    UNSflow.add_celldata!(extra_data, problem.variable_vorticity[1], 
        "NormalVel", map(x->UNSflow.unit(x[1])*x[2], zip(wnormals, wnvel)))
    map(x->UNSflow.add_pointdata!(extra_data, "FlowVel", x[1], x[2]),
        zip(UNSflow.coords(free_wake.geometry), wake_vels))
    mesh = UNSflow.UnstructuredMesh()
    push!(mesh, UNSflow.vorticities(problem))
    UNSflow.add_data!(mesh, extra_data)
    UNSflow.to_vtk_file(mesh, string("output/UVLM_inv_", i))

    influence = deepcopy(UNSflow.vorticities(problem))
    UNSflow.euler!(linked_wake, influence, dt)
    UNSflow.euler!(free_wake, influence, dt)
    UNSflow.increment!(kinem, dt)
    # Transfer geometry...
    linked_wing.geometry.coordinates = map(x->kinem(x), wing_geom.coordinates)
    UNSflow.add_row!(free_wake, 
        linked_wake.geometry.coordinates[1,:], linked_wake.vorticity[1,:])
    linked_wake.geometry.coordinates[2,:] =linked_wake.geometry.coordinates[1,:]
    linked_wake.geometry.coordinates[1,:] = wing.geometry.coordinates[end, :]
end

mesh = UNSflow.UnstructuredMesh()
push!(mesh, problem.variable_vorticity)
println(problem.variable_vorticity)
UNSflow.to_vtk_file(mesh, string("output/UVLM_var_", max+1))
mesh = UNSflow.UnstructuredMesh()
push!(mesh, problem.invariant_vorticity)
println(problem.invariant_vorticity)
UNSflow.to_vtk_file(mesh, string("output/UVLM_inv_", max+1))

end #let

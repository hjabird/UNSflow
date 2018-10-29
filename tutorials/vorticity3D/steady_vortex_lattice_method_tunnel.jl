#===============================================================================
    steady_vortex_lattice_method.jl

    An example of using the vortex lattice method, but this time we also 
    add a tunnel. Variation on steady_vortex_lattice.jl

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
import WriteVTK

let
# STEP 1 : Define a wing. ------------------------------------------------------
# This is done by using equation surf. We take x=-1 as the leading edge, x=1
# as the trailing edge. 
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
spanwise_disc = map(x->cbrt(x), -1:0.1:1)
chordwise_disc = map(x->x, -1 : 0.125 : 1)
wing_geom = UNSflow.discretise(wing_surf, UNSflow.BilinearQuadSurf, 
                                                chordwise_disc, spanwise_disc)

# Lets define some tunnel side walls
offset = 4/3
wallr_fn = x->UNSflow.Vector3D(x[1]*5+3, aspect_ratio * 0.5 * offset, x[2]*2)
walll_fn = x->wallr_fn(x) + UNSflow.Vector3D(0, -offset * aspect_ratio, 0)
wallb_fn = x->UNSflow.Vector3D(x[1]*5+3, aspect_ratio * 0.5 * offset * x[2], -1.7)
wallt_fn = x->UNSflow.Vector3D(x[1]*5+3, aspect_ratio * 0.5 * offset * x[2],  1.7)
wallr_geom = UNSflow.discretise(UNSflow.EquationSurf(wallr_fn),
    UNSflow.BilinearQuadSurf, collect(-1:0.1:1), collect(-1:0.1:1))
walll_geom = UNSflow.discretise(UNSflow.EquationSurf(walll_fn),
    UNSflow.BilinearQuadSurf, collect(-1:0.1:1), collect(-1:0.1:1))
wallb_geom = UNSflow.discretise(UNSflow.EquationSurf(wallb_fn),
    UNSflow.BilinearQuadSurf, collect(-1:0.1:1), collect(-1:0.1:1))
wallt_geom = UNSflow.discretise(UNSflow.EquationSurf(wallt_fn),
    UNSflow.BilinearQuadSurf, collect(-1:0.1:1), collect(-1:0.1:1))

# STEP 2: Define the free stream and wake parameters to setup the problem ------
free_stream = UNSflow.Vector3D(1., 0, 0)
# We have to decide on how to discretise the wake as it travels downstream.
downstream = map(x->x^1.4, collect(0:0.2:3)[2:end] .* 5)
problem = UNSflow.VortexLatticeMethod(wing_geom, free_stream, downstream)
# Add our walls to the problem
var_vorts = problem.problem.variable_vorticities
push!(var_vorts, UNSflow.VortexRingLattice(wallr_geom))
UNSflow.add_bc!(problem.problem, UNSflow.centres(wallr_geom), 
    UNSflow.normals(wallr_geom))
push!(var_vorts, UNSflow.VortexRingLattice(walll_geom))
UNSflow.add_bc!(problem.problem, UNSflow.centres(walll_geom), 
    UNSflow.normals(walll_geom))
push!(var_vorts, UNSflow.VortexRingLattice(wallb_geom))
UNSflow.add_bc!(problem.problem, UNSflow.centres(wallb_geom),
    UNSflow.normals(wallb_geom))
push!(var_vorts, UNSflow.VortexRingLattice(wallt_geom))
UNSflow.add_bc!(problem.problem, UNSflow.centres(wallt_geom), 
    UNSflow.normals(wallt_geom))

# STEP 3: Solve the system -----------------------------------------------------
for i = 1 : 1
    # We call solve to calculate vorticities
    UNSflow.solve!(problem)
    # And then relax_wake! to let the wake curl up as expected.
    UNSflow.relax_wake!(problem; relaxation_factor=0.5)
    # We might want to solve again now the problem geometry has changed...
end

# STEP 4: Post-process ---------------------------------------------------------
# We need the induced velocity to caluclate the forces on the lattice:
ind_vel_fn = x -> 
    UNSflow.induced_velocity(UNSflow.vorticities(problem.problem), x)
# We can calculate some useful stuff using steady_loads:
force, moment, pressure = UNSflow.steady_loads(
    problem.problem.variable_vorticities[1], 
    ind_vel_fn; 
    # The wing vortex lattice is overlapped by that of the wake, so we need to 
    # correct for that - imax indicates the maximum i index of the lattice.
    imax_filament_strs=-problem.problem.variable_vorticities[2].vorticity[1,:])

# Print out some useful things:
println("Wing area is ", UNSflow.area(wing_geom))
println("Forces are ", force)
println("Force coeffs ", 2 * force / 
                        (abs(free_stream)^2 * UNSflow.area(wing_geom)))
println("Moments are ", moment)

# STEP 5: Output to VTK file ---------------------------------------------------
mesh = UNSflow.UnstructuredMesh()
extra_data = UNSflow.MeshDataLinker()
UNSflow.add_celldata!(extra_data, problem.problem.variable_vorticities[1],
    "Pressure", pressure)
push!(mesh, problem.problem.variable_vorticities)
UNSflow.add_data!(mesh, extra_data)
UNSflow.to_vtk_file(mesh, "output/steady_vortex_lattice_tunnel")

end #let

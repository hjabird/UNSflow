
push!(LOAD_PATH, "../../src/")
import UNSflow
import WriteVTK

let

# Make something wing shaped:
surf_fn = x->UNSflow.Vector3D(sqrt(1- x[2]^2) *x[1]/2, x[2]*2, 0)
wing_surf = UNSflow.EquationSurf(surf_fn)
spanwise = map(x->cbrt(x), -1:0.05:1)
wing_geom = UNSflow.discretise(wing_surf, UNSflow.BilinearQuadSurf, collect(-1:0.25:1), spanwise)

# Problem setup
downstream = map(x->x^1.4, collect(0:0.1:3)[2:end] .* 5)
free_stream = UNSflow.Vector3D(1., 0., 0.2)
problem = UNSflow.VortexLatticeMethod(wing_geom, free_stream, downstream)
for i = 1 : 1
    UNSflow.solve!(problem)
    UNSflow.relax_wake!(problem)
end
#=
ind_vel_fn = x -> free_stream + UNSflow.induced_velocity(problem.wake_aerodynamic, x) + UNSflow.induced_velocity(problem.wing_aerodynamic, x)
force, moment, pressure = UNSflow.steady_loads(problem.wing_aerodynamic, ind_vel_fn)


println("Wing area is ", UNSflow.area(wing_geom))
println("Forces are ", force)
println("Moments are ", moment)
=#
# And plot the shape of stuff:
points, cells = UNSflow.to_VtkMesh(convert(Vector{UNSflow.BilinearQuad}, problem.wing_geometry))
points, cells = UNSflow.add_to_VtkMesh(points, cells, convert(Vector{UNSflow.BilinearQuad}, problem.wake_aerodynamic.geometry))
vtkfile = WriteVTK.vtk_grid("output/steady_vortex_lattice", points, cells)
vorticity = vcat(vec(problem.wing_aerodynamic.strengths), vec(problem.wake_aerodynamic.strengths))
WriteVTK.vtk_cell_data(vtkfile, vorticity, "vorticity")
outfiles = WriteVTK.vtk_save(vtkfile)

end #let

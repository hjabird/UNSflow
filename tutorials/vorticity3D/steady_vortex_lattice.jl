
push!(LOAD_PATH, "../../src/")
import UNSflow
import WriteVTK

let

# Make something wing shaped:
surf_fn = x->UNSflow.Vector3D(x[1]/2, x[2]*2, 0)
wing_surf = UNSflow.EquationSurf(surf_fn)
spanwise = vcat(collect(-1:0.02:-0.8), collect(-0.7:0.2:0.7), collect(0.8:0.02:1))
wing_goem = UNSflow.discretise(wing_surf, UNSflow.BilinearQuadSurf, collect(-1:0.2:1), spanwise)

# Problem setup
problem = UNSflow.VortexLatticeMethod(wing_goem, UNSflow.Vector3D(1., 0., 0.2), vcat(collect(0.3:0.3:4), [100]))
for i = 1 : 30
    UNSflow.solve!(problem)
    UNSflow.relax_wake!(problem)
end

# And plot the shape of stuff:
points, cells = UNSflow.to_VtkMesh(convert(Vector{UNSflow.BilinearQuad}, problem.wing_geometry))
points, cells = UNSflow.add_to_VtkMesh(points, cells, convert(Vector{UNSflow.BilinearQuad}, problem.wake_aerodynamic.geometry))
vtkfile = WriteVTK.vtk_grid("output/steady_vortex_lattice", points, cells)
vorticity = vcat(vec(problem.wing_aerodynamic.strengths), vec(problem.wake_aerodynamic.strengths))
WriteVTK.vtk_cell_data(vtkfile, vorticity, "vorticity")
outfiles = WriteVTK.vtk_save(vtkfile)

end #let

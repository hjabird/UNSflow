
import WriteVTK
push!(LOAD_PATH,"../../src/")
import UNSflow

let

a = UNSflow.Lautat3D()
a.strip_divisions = 3
a.strip_discretisation = collect(-1 : 0.05 : 1)
a.strip_centres = collect(-0.5 : 0.5 : 0.5)
a.wing = UNSflow.EquationSurf(x->UNSflow.Vector3D(x[1] * sqrt(1 - 0.9* x[2]^2), x[2] * 4, 0))

UNSflow.discretise_wing_to_geometry!(a)

points, cells = UNSflow.to_VtkMesh(a.wing_geometry[1])
for i = 2 : length(a.wing_geometry)
    points, cells = UNSflow.add_to_VtkMesh(points, cells, a.wing_geometry[i])
end
vtkfile = WriteVTK.vtk_grid("output_LAUTAT3D", points, cells)
outfiles = WriteVTK.vtk_save(vtkfile)



end

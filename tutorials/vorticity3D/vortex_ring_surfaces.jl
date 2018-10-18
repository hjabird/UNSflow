
import WriteVTK
push!(LOAD_PATH,"../../src/")
import UNSflow

let 
wing = UNSflow.EquationSurf(
    x->UNSflow.Vector3D(x[1], x[2], 0))
xdim = 10
ydim = 10
disc_surf = UNSflow.discretise(wing, UNSflow.BilinearQuadSurf, collect(-1:2/xdim:1), collect(-1:2/ydim:1))
ring_lattice = UNSflow.VortexRingLattice(zeros(xdim, ydim), disc_surf)

ext_ring = UNSflow.VortexRing(
    UNSflow.Vector3D(-2,-2,0), UNSflow.Vector3D(2,-2,0), UNSflow.Vector3D(2,2,0), UNSflow.Vector3D(-2,2,0), 1.0
)

surf_coords = vec(UNSflow.centres(disc_surf))
surf_normals = vec(UNSflow.normals(disc_surf))

mat = UNSflow.normal_velocity_influence_matrix(ring_lattice, vec(surf_coords), vec(surf_normals))
known = map(x->UNSflow.dot(UNSflow.induced_velocity(ext_ring, x[1]), x[2]), zip(surf_coords, surf_normals))
solution = mat \ -known
UNSflow.update_using_vorticity_vector!(ring_lattice, solution)


# And output to VTK...
points, cells = UNSflow.to_VtkMesh(convert(Vector{UNSflow.BilinearQuad}, disc_surf))
points, cells = UNSflow.add_to_VtkMesh(points, cells, ext_ring.geometry)
vorticity = vec(ring_lattice.strengths)
push!(vorticity, 1.0) #For our vortex ring.

vtkfile = WriteVTK.vtk_grid("output/lattice_in_ring", points, cells)
WriteVTK.vtk_cell_data(vtkfile, vorticity, "vorticity")
outfiles = WriteVTK.vtk_save(vtkfile)

end # End let
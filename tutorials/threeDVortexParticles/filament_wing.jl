#===============================================================================
orbiting_a_vortex_filament.jl

A few particles orbiting round a vortex filament.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow.
using Dierckx
using SpecialFunctions
import Base.size
include("../../src/kinem.jl")
include("../../src/utils.jl")
include("../../src/delVort.jl")
include("../../src/lowOrder2D/typedefs.jl")
include("../../src/lowOrder2D/calcs.jl")
include("../../src/lowOrder3D/typedefs.jl")
include("../../src/lowOrder3D/calcs.jl")
include("../../src/lowOrder3D/3DLautat.jl")
using WriteVTK  # We'll use this package to output to VTK for visualisation.

y = [-2., -1., 0., 1., 2.]
c = [0.7, 0.9, 1., 0.9, 0.7]
n_chords = 5
chords = Vector{WingChordSection}(n_chords)
for i = 1 : n_chords
    chords[i] = WingChordSection()
    chords[i].camber_line = Spline1D([-1., 0.0, 1.], [0.0, 0.1, 0.0], k=2)
    chords[i].LE_location.y = y[i]
    chords[i].TE_location.y = y[i]
    chords[i].LE_location.x = 0.5 - c[i] / 2.0
    chords[i].TE_location.x = 0.5 + c[i] / 2.0
end # Good.

wing = StripDefinedWing(
    chords,
    ThreeDVector(0.25, 2.5, 0),
    ThreeDVector(0.75, 2.5, 0),
    ThreeDVector(0.25, -2.5, 0),
    ThreeDVector(0.75, -2.5, 0),
)

filament_positions = Vector{Float64}([-cos(x) for x in 0 : 0.1 : pi])
fil_wing = build_vortex_filament_wing_geometry(
    wing, filament_positions
)
add_vorticity!(wing, filament_positions, fil_wing, 2, x->x)
for i = 1 : length(wing.strips)
    add_vorticity!(wing, filament_positions, fil_wing, i, x->cos(i * pi * x))
end
k_sloc = kelvin_particles_span_shedding_locations(wing, 0.1)
k_sind = k_particle_shedding_locs_to_bv_index(k_sloc, length(wing.strips))
free_stream = x->ThreeDVector(1.0, 0, 0)
particles = shed_initial_particles(wing, free_stream, k_sloc, 0.1, 0.1)
vf = get_particle_vorticity_function(particles, k_sind, 0.1, x->ThreeDVector(1, 0, 0))
a = vf(wing_to_bv_vector(fil_wing), zeros(10))
for i = 1 : length(a)
    particles[i].vorticity *= a[i]
end
particles = vcat(particles, shed_particles(wing, k_sloc, 0.1, particles))
sim = VortexParticleWakeLAUTATSolution()
sim.wing = wing
sim.filament_wing = fil_wing



fils = convert(Vector{ThreeDStraightVortexFilament}, fil_wing)
nfils = length(fils)

#points = zeros(3, nfils * 2)
#cells = Array{MeshCell, 1}(nfils)
#celldata = zeros(nfils)
#for j = 1 : nfils
    #points[:, 2 * j] = convert(Array{Float64, 1}, fils[j].start_coord)
    #points[:, 2 * j - 1] = convert(Array{Float64, 1}, fils[j].end_coord)
    #cells[j] = MeshCell(VTKCellTypes.VTK_LINE, [2 * j, 2 * j - 1])
    #celldata[j] = fils[j].vorticity
#end
#vtkfile = vtk_grid(string("output/Hello_filaments", 0), points, cells)
#vtk_cell_data(vtkfile, celldata, "vorticity")
#outfiles = vtk_save(vtkfile)

npart = length(particles)
points2 = zeros(3, npart)
cells2 = Array{MeshCell, 1}(npart)
vort2 = zeros(3, npart)
for j = 1 : npart
    points2[:, j] = convert(Array{Float64, 1}, particles[j].coord)
    cells2[j] = MeshCell(VTKCellTypes.VTK_VERTEX, [j])
    vort2[:, j] = convert(Array{Float64, 1}, particles[j].vorticity)
end
vtkfile2 = vtk_grid(string("output/Hello_filaments2", 0), points2, cells2)
vtk_cell_data(vtkfile2, vort2, "vorticity")
outfiles = vtk_save(vtkfile2)

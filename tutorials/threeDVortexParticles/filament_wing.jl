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
import SpecialFunctions
import Base.size
include("../../src/kinem.jl")
include("../../src/utils.jl")
include("../../src/delVort.jl")
include("../../src/lowOrder2D/typedefs.jl")
include("../../src/lowOrder2D/calcs.jl")
include("../../src/lowOrder3D/typedefs.jl")
include("../../src/lowOrder3D/calcs.jl")
include("../../src/lowOrder3D/3DLautat.jl")
import WriteVTK  # We'll use this package to output to VTK for visualisation.

y = [-2., 0., 2.]
c = [0.7, 1., 0.7]
n_chords = 3
chords = Vector{WingChordSection}(n_chords)
for i = 1 : n_chords
    chords[i] = WingChordSection()
    chords[i].camber_line = Spline1D([-1., 0.0, 1.], [0.0, 0.0, 0.0], k=2)
    chords[i].LE_location.y = y[i]
    chords[i].TE_location.y = y[i]
    chords[i].LE_location.x = 0.5 - c[i] / 2.0
    chords[i].TE_location.x = 0.5 + c[i] / 2.0
    chords[i].LE_location.z = y[i]^2 * 0.1
    chords[i].TE_location.z = y[i]^2 * 0.1
end # Good.

t0wing = StripDefinedWing(
    chords,
    ThreeDVector(0.25, 2.5, 0.8),
    ThreeDVector(0.75, 2.5, 0.8),
    ThreeDVector(0.25, -2.5, 0.8),
    ThreeDVector(0.75, -2.5, 0.8),
)

filament_positions = Vector{Float64}([-cos(x) for x in 0 : 0.1 : pi])
original_fil_wing = build_vortex_filament_wing_geometry(
    t0wing, filament_positions
    )
fil_wing = deepcopy(original_fil_wing)

kinem = ThreeDCoordinateTransform((x,t)->x + sin(t) * ThreeDVector(0,0,0))
k_sloc = kelvin_particles_span_shedding_locations(t0wing, 0.1)
k_sind = k_particle_shedding_locs_to_bv_index(k_sloc, length(t0wing.strips))
free_stream = x->ThreeDVector(1.0, 0, 0.0)
wake = ThreeDVortexParticleSet()
n_fourier_terms = 1
dt = 0.1
old_bvs = zeros(length(t0wing.strips) * 2)

fil_wing = transform_ThreeDVectors(x->kinem(x), original_fil_wing)
particles = shed_initial_particles(t0wing, free_stream, k_sloc, dt, 0.1)
ind_vel_external = x->free_stream(x)
f_terms = solve_new_vortex_particle_vorticities_and_assign!(
    t0wing, fil_wing, kinem,
    filament_positions, old_bvs, abs(free_stream(ThreeDVector(0,0,0))),
    particles, k_sind, dt, ind_vel_external, n_fourier_terms
    )
set_wing_to_fourier_set!(fil_wing, t0wing, filament_positions, f_terms)
old_bvs = wing_to_bv_vector(fil_wing)
wake += particles

for n = 1 : 2
    # Convection
    e_ind_vel = x->free_stream(x) + ind_vel(fil_wing, x)
    e_ind_dvort = x->ThreeDVector(0,0,0) #ind_dvortdt(x, fil_wing)
    euler_forward_step!(wake, e_ind_vel, e_ind_dvort, dt)
    increment!(kinem, dt)
    # Shedding and variable definition updates
    fil_wing = transform_ThreeDVectors(x->kinem(x), original_fil_wing)
    wing = transform_ThreeDVectors(x->kinem(x), t0wing)

    last_shed_particles = deepcopy(particles)
    old_bvs = wing_to_bv_vector(fil_wing)
    particles = shed_particles(wing, k_sloc, 0.1, last_shed_particles)
    ind_vel_external = x->ind_vel(wake, x) + free_stream(x)
    # Finding the new vorticities
    f_terms = solve_new_vortex_particle_vorticities_and_assign!(
        t0wing, fil_wing, kinem,
        filament_positions, old_bvs, abs(free_stream(ThreeDVector(0,0,0))),
        particles, k_sind, dt, ind_vel_external, n_fourier_terms
        )
    wake += particles
    print(string("Iteration ", n, ". ", length(wake), " vortex particles.\n"))
end

# VKT EXPORT ===================================================================
fils = convert(Vector{ThreeDStraightVortexFilament}, fil_wing)
nfils = length(fils)

points = zeros(3, nfils * 2)
cells = Array{WriteVTK.MeshCell, 1}(nfils)
celldata = zeros(3, nfils)
for j = 1 : nfils
    points[:, 2 * j] = convert(Array{Float64, 1}, fils[j].start_coord)
    points[:, 2 * j - 1] = convert(Array{Float64, 1}, fils[j].end_coord)
    cells[j] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE, [2 * j, 2 * j - 1])
    celldata[:, j] = convert(Array{Float64, 1}, ThreeDVector(fils[j].vorticity, 0, 0))
end
vtkfile = WriteVTK.vtk_grid(string("output/Hello_filaments", 0), points, cells)
WriteVTK.vtk_cell_data(vtkfile, celldata, "vorticity")
outfiles = WriteVTK.vtk_save(vtkfile)

npart = length(wake)
points2 = zeros(3, npart)
cells2 = Array{WriteVTK.MeshCell, 1}(npart)
vort2 = zeros(3, npart)
for j = 1 : npart
    points2[:, j] = convert(Array{Float64, 1}, wake[j].coord)
    cells2[j] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [j])
    vort2[:, j] = convert(Array{Float64, 1}, wake[j].vorticity)
end
vtkfile2 = WriteVTK.vtk_grid(string("output/Hello_filaments2", 0), points2, cells2)
WriteVTK.vtk_cell_data(vtkfile2, vort2, "vorticity")
outfiles = WriteVTK.vtk_save(vtkfile2)

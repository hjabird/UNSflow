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

y_span = 1
y_values = y_span * sin.(linspace(-pi/2, pi/2, 7))
print(string("y_values: ", y_values, "\n"))
c_fn = x -> 1 #- 0.2* x ^2
c_values = c_fn.(y_values / y_span)
z_fn = z -> 0
z_values = z_fn.(y_values / y_span)
chords = Vector{WingChordSection}(length(y_values) - 2)
for i = 2 : length(y_values) - 1
    chords[i-1] = WingChordSection()
    chords[i-1].camber_line = Spline1D([-1., 0.0, 1.], [0.0, 0.0, 0.0], k=2)
    chords[i-1].LE_location.y = y_values[i]
    chords[i-1].TE_location.y = y_values[i]
    chords[i-1].LE_location.x = 0.5 - c_values[i] / 2.0
    chords[i-1].TE_location.x = 0.5 + c_values[i] / 2.0
    chords[i-1].LE_location.z = z_values[i]
    chords[i-1].TE_location.z = z_values[i]
end # Good.

t0wing = StripDefinedWing(
    chords,
    ThreeDVector(0.5 - c_values[end] / 2.0, y_values[end], z_values[end]),
    ThreeDVector(0.5 + c_values[end] / 2.0, y_values[end], z_values[end]),
    ThreeDVector(0.5 - c_values[1] / 2.0, y_values[1], z_values[1]),
    ThreeDVector(0.5 + c_values[1] / 2.0, y_values[1], z_values[1]),
)

#filament_positions = Vector{Float64}(-1 : 0.01 : 1)[2:end-1]
filament_positions = map(x->-cos(x), Vector{Float64}(-0 : 0.1 : pi)[2:end])
original_fil_wing = build_vortex_filament_wing_geometry(
    t0wing, filament_positions
    )
fil_wing = deepcopy(original_fil_wing)

kinem = ThreeDCoordinateTransform((x,t)->rotate_about_y(x, 0.4))
k_sloc = kelvin_particles_span_shedding_locations(t0wing, 0.1)
k_sind = k_particle_shedding_locs_to_bv_index(k_sloc, length(t0wing.strips))
free_stream = x->ThreeDVector(0.921061, 0.0, 0.0) #0.4 rad
wake = ThreeDVortexParticleSet()
n_fourier_terms = 2
dt = 0.025
old_bvs = zeros(length(t0wing.strips) * 2)

fil_wing = transform_ThreeDVectors(x->kinem(x), original_fil_wing)
old_fil_wing = fil_wing
fil_wing = transform_ThreeDVectors(x->kinem(x), old_fil_wing)
wing = transform_ThreeDVectors(x->kinem(x), t0wing)
particles = shed_initial_particles(wing, free_stream, k_sloc, dt, 0.1)
ind_vel_external = x->free_stream(x)
f_terms = solve_new_vortex_particle_vorticities_and_assign!(
    t0wing, fil_wing, kinem,
    filament_positions, old_bvs, abs(free_stream(ThreeDVector(0,0,0))),
    particles, k_sind, dt, ind_vel_external, n_fourier_terms
    )
set_wing_to_fourier_set!(fil_wing, t0wing, filament_positions, f_terms)
old_bvs = wing_to_bv_vector(fil_wing)
wake += particles
print(string("Fourier terms: ", f_terms, "\n"))

for n = 1 : 90
    # Convection
    e_ind_vel = x->free_stream(x) + ind_vel(fil_wing, x)
    e_ind_dvort = x->ThreeDVector(0,0,0) #ind_dvortdt(x, fil_wing)
    euler_forward_step!(wake, e_ind_vel, e_ind_dvort, dt)
    increment!(kinem, dt)
    # Shedding and variable definition updates
    old_fil_wing = fil_wing
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
    print(string("Fourier terms: ", f_terms, "\n"))
end

fil_wing = transform_ThreeDVectors(x->kinem(x), original_fil_wing)
ind_vel_p = x->ind_vel_external(x) + ind_vel(fil_wing, x) + ind_vel(wake, x)
pd_fn = pressure_distribution(fil_wing, old_fil_wing, t0wing, filament_positions, dt,
    ind_vel_p, f_terms, 1.0)
pd_ys = [j for j = 0.75:0.5:length(t0wing.strips)+0.25, i = 1 : length(filament_positions)]
pd_xs = [filament_positions[i] for j = 0.75:0.5:length(t0wing.strips)+0.25, i = 1 : length(filament_positions)]
pdist = pd_fn.(pd_xs, pd_ys)
print(pdist)

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

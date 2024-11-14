#### Plotting the results is not discussed here. For visualizing the results, you need to go to the GUI version of Tidy3D ##########

!pip install tidy3d
import tidy3d.web as web
web.configure("wclpAIGF4JlvJlABZs9BOY3yg53fsAAd54VzLiVcCFPu7pFC")

!pip install gdstk

import tidy3d.web as web # if needed
web.test()

import numpy as np
import matplotlib.pyplot as plt
import tidy3d as td
import tidy3d.web as web
from tidy3d.plugins.resonance import ResonanceFinder

import gdstk

lambda_min = 0.8
lambda_max = 1.1

monitor_lambda = np.linspace(lambda_min, lambda_max, 101)
monitor_freq = td.constants.C_0/monitor_lambda

dl = 0.02 # uniform sampling mesh size

sim_size = Lx, Ly, Lz = (3.5, 3.5, 2)
offset_monitor = 0.4 # monitor position in the vertical direction

t_start = 0.2e-12  # from inspection
t_stop = 4e-12

vacuum = td.Medium(permittivity = 1, name = 'vacuum')
GaAs_permittivity = 3.55**2
GaAs = td.Medium(permittivity = GaAs_permittivity, name = 'GaAs')

freq0 = 302.675e12
fwidth = 144.131e12

slab_side_length = 10
slab_height = 0.2

slab = td.Structure(
    geometry = td.Box(
        center = (0, 0, 0),
        size = (slab_side_length, slab_side_length, slab_height)
    ),
    medium = GaAs,
    name = 'GaAs slab'
)

lib = gdstk.Library()

# Geometry must be placed in cells.
cell = lib.new_cell("BULLSEYE")

central_inner_radius = 0.355                   #initial inner radius
ring_width = 0.07                     #rings' width
ring_height = 0.2                   #rings' height
theta1 = 3 * np.pi / 180
theta2 = 87 * np.pi / 180
ringsNumbers = 7
vertexIncrement = 0.1075
ringNumber = 7
xOrigin = 0
yOrigin = 0


for i in range(ringsNumbers):
    for j in range(4):
        gcSlice = gdstk.ellipse(
            (xOrigin, yOrigin),
            central_inner_radius + i * (vertexIncrement + ring_width),
            central_inner_radius + i * (vertexIncrement + ring_width) + ring_width,
            theta1 + j*np.pi/2,
            theta2 + j*np.pi/2,
            layer=1,
            datatype=1,
            tolerance=0.0005
        )
        cell.add(gcSlice)
gc_etch = td.Geometry.from_gds(
        cell, gds_layer=1, axis=2, slab_bounds=(-ring_height/2, ring_height/2)
    )

gc_etch.plot(z=0)
plt.show()

mat_etch = td.Medium(permittivity = 1, name = 'air')
etch = td.Structure(
    geometry = gc_etch,
    medium = mat_etch,
    name = 'etch'
)

dipole_source = td.PointDipole(
    center = (0, 0, 0),
    source_time = td.GaussianPulse(
        freq0 = freq0,
        fwidth = fwidth),
    polarization = 'Ex',
    name = 'dipole_source'
)


far_field_monitor = td.FieldProjectionAngleMonitor(
    center = (0, 0, offset_monitor),
    size = (td.inf, td.inf, 0),
    name = 'far_field_monitor',
    freqs = monitor_freq,
    normal_dir = '+',
    phi = np.linspace(0, 2 * np.pi, 181),
    theta = np.linspace(0, np.pi, 91)
)

sim_far = td.Simulation(
    center = (0, 0, 0),
    size = sim_size,
    grid_spec = td.GridSpec.uniform(dl = dl),
    run_time = t_stop,
    sources = [dipole_source],
    monitors = [far_field_monitor],
    structures = [slab, etch],
    medium = vacuum,
    shutoff = 1e-7,
    symmetry = [-1, 1, 0]
)

sim_far.plot_3d()

task_id = web.upload(sim_far, task_name = 'cavity_far')

web.start(task_id)
web.monitor(task_id, verbose = True)

sim_data = web.load(task_id, path='data/sim_data_far.hdf5')
print(sim_data.log)

sim_data = td.SimulationData.from_file(fname = 'data/sim_data_far.hdf5')

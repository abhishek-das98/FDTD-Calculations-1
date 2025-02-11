


!pip install tidy3d
import tidy3d.web as web
web.configure("wclpAIGF4JlvJlABZs9BOY3yg53fsAAd54VzLiVcCFPu7pFC") # API Key for my personal account on Tidy3D, This will vary for different acccounts

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

t_start = 0.2e-12  # from inspection, you can check for different t_start if the quality factor changes at all
t_stop = 4e-12

vacuum = td.Medium(permittivity = 1, name = 'vacuum') # For inverse design calculations, take the permittivity to be a little higher
GaAs_permittivity = 3.55**2
GaAs = td.Medium(permittivity = GaAs_permittivity, name = 'GaAs')

freq0 = 302.675e12 # taken from the Lumerical file supplied by Dr. Dima
fwidth = 144.131e12 # taken from the Lumerical file supplued by Dr. Dima

slab_side_length = 10 # arbitrary, ultimately cuts down to the simulation size
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

central_inner_radius = 0.355           #initial inner radius
ring_width = 0.07                    #rings' width
ring_height = 0.2                   #rings' height
theta1 = 3 * np.pi / 180            # bridge angle is 3 degrees
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

time_monitor_1 = td.FieldTimeMonitor(
    name = 'time_monitor_1',
    size = [0, 0, 0],
    center = [0.1, 0.1, 0],
    start = t_start
)

time_monitor_2 = td.FieldTimeMonitor(
    name = 'time_monitor_2',
    size = [0, 0, 0],
    center = [0.1, 0.2, 0],
    start = t_start
)

time_monitor_3 = td.FieldTimeMonitor(
    name = 'time_monitor_3',
    size = [0, 0, 0],
    center = [0.2, 0.2, 0],
    start = t_start
)

time_monitor_4 = td.FieldTimeMonitor(
    name = 'time_monitor_4',
    size = [0, 0, 0],
    center = [0.2, 0.1, 0],
    start = t_start
)

time_monitors = [time_monitor_1, time_monitor_2, time_monitor_3, time_monitor_4]

sim_resonance = td.Simulation(
    center = (0, 0, 0),
    size = sim_size,
    grid_spec = td.GridSpec.uniform(dl = dl),
    run_time = t_stop,
    sources = [dipole_source],
    monitors = time_monitors,
    structures = [slab, etch],
    medium = vacuum,
    shutoff = 1e-7,
    symmetry = [-1, 1, 0]   # PEC Symmetry along the x-direction, PMC symmetry along the y-direction; This will be applicable for x-polarized light
)

sim_resonance.plot_3d()

task_id = web.upload(sim_resonance, task_name = 'cavity_resonance_target_1000')

web.start(task_id)
web.monitor(task_id, verbose = True)

sim_data = web.load(task_id, path='data/sim_data_resonance_target_1000.hdf5')
print(sim_data.log)

sim_data = td.SimulationData.from_file(fname = 'data/sim_data_resonance_target_1000.hdf5')

resonance_finder = ResonanceFinder(freq_window=(monitor_freq[-1], monitor_freq[0]))
resonance_data = resonance_finder.run(signals=sim_data.data)
resonance_data.to_dataframe()

import numpy as np
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(1, 1, tight_layout=True, figsize=(7, 4))

# Extract the time response data
time_response = sim_data["time_monitor_1"].Ex.squeeze()

# Perform FFT and get the frequency response
freq_response = np.abs(np.fft.fft(time_response))

# Calculate the frequency values
freqs = np.linspace(0, 1 / sim_data.simulation.dt, len(time_response))

# Calculate the corresponding wavelengths
wavelengths = td.constants.C_0 / freqs

# Select the indices within the monitor frequency range
plot_inds = np.where((monitor_freq[-1] < freqs) & (freqs < monitor_freq[0]))

# Plot the frequency response with respect to wavelength
ax1.plot(wavelengths[plot_inds], freq_response[plot_inds])
ax1.set_xlabel("Wavelength (m)")
ax1.set_ylabel("Amplitude")
ax1.set_title("Frequency Response vs. Wavelength")
ax1.grid(True)

# Reverse the x-axis to show decreasing wavelength (optional)
#ax1.set_xlim(ax1.get_xlim()[::-1])

plt.show()
plt.savefig('Resonance_data')

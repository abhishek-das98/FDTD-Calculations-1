pip install legume-gme

import legume

import numpy as np
import matplotlib.pyplot as  plt
import time

from legume.minimize import Minimize

%load_ext autoreload
%autoreload 2

central_radius = 0.3550
slab_height = 0.2
ring_width = 0.07
ring_number = 7
ring_period = 0.155

GaAs_permittivity = 12.25

Lx, Ly = 6, 6

def rings(dcenter, dperiod): # The function calculates the structure given the deviations from the initial central ring radius and ring period


  lattice = legume.Lattice([Lx, 0], [0, Ly])
  phc = legume.PhotCryst(lattice)

  phc.add_layer(d = slab_height, eps_b = GaAs_permittivity)

  angles = np.linspace(0, 2*np.pi, 361)

  rings = []
  for j in range(ring_number):

    inner_radius = central_radius + dcenter + (ring_width + ring_period + dperiod) * j
    outer_radius = inner_radius + ring_width
    xs = []
    ys = []
    for angle in angles:
      xs.append(outer_radius * np.cos(angle))
      ys.append(outer_radius * np.sin(angle))
    for angleR in angles[::-1]:
      xs.append(inner_radius * np.cos(angleR))
      ys.append(inner_radius * np.sin(angleR))
    ring = legume.Poly(
        eps = 1.0,
        x_edges = xs,
        y_edges = ys
    )
    rings.append(ring)

  phc.add_shape(rings)

  return phc

# For visualizing the structure
phc = rings(0.0, 0.0)
legume.viz.structure(phc, yz=True, figsize=4., cbar=True)

# This function calculates the kpoints given the number of k points and lattice constants
def get_kpoints(Lx, Ly, nkx, nky):

  # sample nkx and nky points in {kx, ky} space in a uniform grid
  # spacing between two reciprocal vectors is 2*pi/N

  kx = np.linspace(0, (nkx-1)/nkx*2*np.pi/Lx, nkx)
  ky = np.linspace(0, (nky-1)/nky*2*np.pi/Ly, nky)
  kxg, kyg = np.meshgrid(kx, ky)
  kxg = kxg.ravel()
  kyg = kyg.ravel()

  kpoints = np.vstack((kxg, kyg))

  return kpoints


def gme_cavity(dx, dy, gmax, truncate_g, options): # function for running the gme simulation

  bullseye = rings(dx, dy)

  options['compute_im'] = False      
 
  gme = legume.GuidedModeExp(bullseye, gmax = gmax, truncate_g = truncate_g)

  kpoints_number = 2
  kpoints = get_kpoints(Lx, Ly, kpoints_number, kpoints_number)

  gme.run(kpoints = kpoints, **options)

  (freq_im, _, _) = gme.compute_rad(0, [1100])

  Q = gme.freqs[0, 1100]/2/freq_im[0]

  return(gme, Q)

options = {
    'gmode_inds': [0],
    'verbose': True,
    'numeig': 1500,
    #'eig_sigma': 1.07,
    'gradients': 'approx'
}

gmax = 4

truncate_g = 'abs'
dx = 0.0
dy = 0.0

(gme, Q) = gme_cavity(dx, dy, gmax, truncate_g, options)

# Print the computed quality factor
print("Cavity quality factor: %1.2f" %Q)

# We can also visualize the cavity and the mode profile of the fundamental mode

ax = legume.viz.field(gme, 'e', 0, 1100, z=slab_height/2, component='x', val='abs', N1=300, N2=200)




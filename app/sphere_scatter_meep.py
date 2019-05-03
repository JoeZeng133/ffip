# %% modules import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ffip
import json
import meep as mp
from meep.materials import Au
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def get_face_volumes(vol=mp.Volume(mp.Vector3())):
    return [
        mp.Volume(vol.center + mp.Vector3(x=-vol.size.x/2),
                  size=mp.Vector3(y=vol.size.y, z=vol.size.z)),
        mp.Volume(vol.center + mp.Vector3(x=+vol.size.x/2),
                  size=mp.Vector3(y=vol.size.y, z=vol.size.z)),
        mp.Volume(vol.center + mp.Vector3(y=-vol.size.y/2),
                  size=mp.Vector3(x=vol.size.x, z=vol.size.z)),
        mp.Volume(vol.center + mp.Vector3(y=+vol.size.y/2),
                  size=mp.Vector3(x=vol.size.x, z=vol.size.z)),
        mp.Volume(vol.center + mp.Vector3(z=-vol.size.z/2),
                  size=mp.Vector3(y=vol.size.y, x=vol.size.x)),
        mp.Volume(vol.center + mp.Vector3(z=+vol.size.z/2),
                  size=mp.Vector3(y=vol.size.y, x=vol.size.x))
    ]


# %%
dx = 2.0e-3
dt = dx * 0.5
resolution = 1/dx
fmin = 1/1000.0e-3
fmax = 1/300.0e-3
fcen = (fmin + fmax) / 2
df = fmax - fmin
dpml = 15.0e-3
nfreq = 20
r = 30e-3

cell_size = mp.Vector3(100.0e-3, 100.0e-3, 100.0e-3)
pml_layers = [mp.PML(thickness=dpml)]
geometry = [mp.Sphere(radius=r, material=Au)]
src = [mp.Source(
    mp.GaussianSource(fcen, fwidth=df),
    component=mp.Ex,
    center=mp.Vector3(0, 0, -cell_size.z/2+dpml),
    size=mp.Vector3(cell_size.x, cell_size.y, 0)
)
]

# %%
sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    sources=src,
    resolution=resolution,
    k_point=mp.Vector3()
)

flux_box_size = mp.Vector3(65e-3, 65e-3, 65e-3)
flux_volumes = get_face_volumes(
    mp.Volume(center=mp.Vector3(), size=flux_box_size))
flux_regions = []
for i in range(6):
    flux_regions.append(
        mp.FluxRegion(
            volume=flux_volumes[i], weight=-1 if i % 2 == 0 else 1
        )
    )

inc_face = sim.add_flux(fcen, df, nfreq, flux_regions[4])
inc_box = sim.add_flux(fcen, df, nfreq, *flux_regions)

# %% run and save files on hard-drive
c = input("to skip running simulation, type 1")
if c != '1':
    sim.run(until=dt*1000)
    sim.save_flux('inc_face', inc_face)
    sim.save_flux('inc_box', inc_box)

# %% retrive data from files
sim.load_flux('inc_face', inc_face)
sim.load_flux('inc_box', inc_box)
inc_flux = np.asarray(mp.get_fluxes(inc_face))

# %% configure the second simulation
sim.reset_meep()
sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    sources=src,
    resolution=resolution,
    geometry=geometry,
    k_point=mp.Vector3()
)

tot_box = sim.add_flux(fcen, df, nfreq, *flux_regions)

# %% run the second simulation
c = input("to skip running simulation, type 1")
if c != '1':
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            c=mp.Ex,
            pt=mp.Vector3(z=-flux_box_size.z/2),
            decay_by=1e-4,
            dt=1/(fcen * 2 * pi)
        )
    )
    # sim.run(until=dt*10000)
    sim.save_flux('tot_box', tot_box)

# %% load data from the second simulation
sim.load_flux('tot_box', tot_box)
tot_flux = np.asarray(mp.get_fluxes(tot_box))
freqs = np.asarray(mp.get_flux_freqs(tot_box))

# %%
plt.plot(1/freqs, tot_flux / inc_flux)
plt.title('absorption vs wavelength')

plt.show()

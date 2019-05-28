#%% modules import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# %load_ext autoreload
# %autoreload 1
# %aimport ffip

#%%
dx = 43.2e-3/ffip.um_scale
sc = 0.5
dt = sc * dx

lmin = 500e-3/ffip.um_scale
lmax = 2000e-3/ffip.um_scale
ls = (lmin + lmax) / 2
l = np.linspace(lmax, lmin, 100)
ft = 1 / l
fs = 1 / ls
d = 400e-3/ffip.um_scale

omega = 2 * pi * ft
material = ffip.FePt(fs)
er = material.get_epsilon(ft)
ur = 1
eta = np.sqrt(ur / er)
k = omega * np.sqrt(ur * er)

R1 = (eta - 1) / (eta + 1)
T1 = R1 + 1
R2 = (1 - eta) / (eta + 1)
T2 = R2 + 1

# total reflection and transmission
Rt = (R1 + R2 * np.exp(-2j * k * d)) / (1 + R1 * R2 * np.exp(-2j * k * d))
Tt = (T1 * T2 * np.exp(-1j * k * d)) / (1 + R1 * R2 * np.exp(-2j * k * d))

size = ffip.Vector3(10, 15, 42) * dx

src_func = ffip.Gaussian1(frequency=fs)

pmls = [ffip.PML(thickness=6*dx,direction='z')]

geometry = [ffip.Box(size=ffip.Vector3(1e9, 1e9, d), material=material)]

sources = [ffip.Source(
    function=src_func, 
    center=ffip.Vector3(0, 0, -15*dx), 
    size=ffip.Vector3(size.x, size.y, 0), 
    amplitude=1,
    field_component='Ex') ]

periodic = [1, 1, 0]

output_filename1 = 'output1.h5'
output_filename2 = 'output2.h5'

sim1 = ffip.Simulation(
    size=size, 
    resolution=1/dx,
    pmls=pmls,
    sources=sources,
    fields_output_file=output_filename1,
    periodic=periodic
)

sim2 = ffip.Simulation(
    size=size,
    resolution=1/dx,
    pmls=pmls,
    sources=sources,
    geometry=geometry,
    fields_output_file=output_filename2,
    periodic=periodic
)

inc_dft = sim1.add_flux_region(
    center=ffip.Vector3(0, 0, -10*dx),
    size=ffip.Vector3(5*dx, 5*dx, 0), 
    frequency=ft
)

tot_dft = sim2.add_flux_region(
    center=ffip.Vector3(0, 0, -10*dx),
    size=ffip.Vector3(5*dx, 5*dx, 0), 
    frequency=ft
)

trn_dft = sim2.add_flux_region(
    center=ffip.Vector3(0, 0, 10*dx),
    size=ffip.Vector3(5*dx, 5*dx, 0), 
    frequency=ft
)

#%%
sim1.run(stop_condition=ffip.run_until_time(src_func.end_time*10), np=2)
sim2.run(stop_condition=ffip.run_until_time(src_func.end_time*10), np=2)

#%%
pts = np.stack(np.meshgrid(ft, 0, 0, 0, indexing='ij'), axis=-1)

inc = inc_dft.get_values(scale=2)
ref = tot_dft.get_values(substract=inc_dft, scale=2)
trn = trn_dft.get_values(scale=2)

R = np.sqrt(-ref / inc)
T = np.sqrt(trn / inc)

#%%
plt.figure(1)
plt.plot(l, np.abs(Rt), linewidth=0.5)
plt.plot(l, np.abs(R), 's', markersize=1)
plt.title('Reflectance vs Wavelength')

plt.figure(2)
plt.plot(l, np.abs(Tt), linewidth=0.5)
plt.plot(l, np.abs(T), 's', markersize=1)
plt.title('Transmittance vs Wavelength')

plt.show()

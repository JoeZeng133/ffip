#%% modules import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#%%
dx = 2e-3/ffip.um_scale
sc = 0.5
dt = dx * sc
lmin = 200e-3/ffip.um_scale
lmax = 1000e-3/ffip.um_scale
fcen = (1/lmin + 1/lmax)/2
d = 40e-3/ffip.um_scale
dpml = 10e-3/ffip.um_scale

l = np.linspace(lmax, lmin, 100)
ft = 1 / l
omega = 2 * pi * ft
material = ffip.Au
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

size = ffip.Vector3(12e-3/ffip.um_scale, 12e-3/ffip.um_scale, 80e-3/ffip.um_scale)

src_func = ffip.Gaussian1(frequency=fcen)

pmls = [ffip.PML(thickness=dpml,direction='z')]

geometry = [ffip.Block(size=ffip.Vector3(1e9, 1e9, d), material=material)]

sources = [ffip.Source(
    function=src_func, 
    center=ffip.Vector3(0, 0, -size.z/2+dpml), 
    size=ffip.Vector3(size.x, size.y, 0), 
    amplitude=1/dx,
    field_component='Ex') ]

periodic = [1, 1, 0]

output_filename1 = 'output3.h5'
output_filename2 = 'output4.h5'

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

mon_region = sim1.add_dft_fields(size=ffip.Vector3(size.x, size.y, size.z), frequency=[fcen], field_component='Ex')
interval = src_func.end_time
pt = ffip.Vector3(0, 0, size.z/2-dpml)
# sim1.run(stop_condition=ffip.run_until_time(interval*4))
sim1.run(stop_condition=ffip.run_until_fields_decay(position=pt,field_component='Ex', time_interval_examined=interval))

# interp = mon_region.get_interpolant(method='nearest', fill_value=None)
# fp = fcen
# xp = np.linspace(-size.x/2, size.x/2, 80)
# yp = np.linspace(-size.y/2, size.y/2, 80)
# zp = np.linspace(-size.z/2, size.z/2, 80)
# pts = np.stack(np.meshgrid(fp, zp, yp, xp, indexing='ij'), axis=-1)
# val = np.squeeze(interp(pts))


input("debug ended")

#%%
inc_dft = sim1.add_flux_region(
    center=ffip.Vector3(0, 0, -25e-3/ffip.um_scale),
    size=ffip.Vector3(5e-3/ffip.um_scale, 5e-3/ffip.um_scale, 0), 
    frequency=ft
)

tot_dft = sim2.add_flux_region(
    center=ffip.Vector3(0, 0, -25e-3/ffip.um_scale),
    size=ffip.Vector3(5e-3/ffip.um_scale, 5e-3/ffip.um_scale, 0), 
    frequency=ft
)

trn_dft = sim2.add_flux_region(
    center=ffip.Vector3(0, 0, 25e-3/ffip.um_scale),
    size=ffip.Vector3(5e-3/ffip.um_scale, 5e-3/ffip.um_scale, 0), 
    frequency=ft
)

#%%
interval = src_func.end_time / 2
sim1.run(stop_condition=ffip.run_until_time(time=1e3*dt))
sim2.run(stop_condition=ffip.run_until_fields_decay(field_component='Ex', time_interval_examined=interval))
# sim2.run(stop_condition=ffip.run_until_time(1e4*dt))

# #%%
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

#%% configure simulation
import numpy as np
import matplotlib.pyplot as plt
import h5py
from .. import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

dx = 43.2
sc = 0.5
dt = sc * dx

lmin = 500
lmax = 2000
ls = (lmin + lmax) / 2
l = np.linspace(lmax, lmin, 100)
ft = 1 / l
fs = 1 / ls
d = 400

omega = 2 * pi * ft
material = ffip.fic
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

geometry = [ffip.Block(size=ffip.Vector3(1e9, 1e9, d), material=material)]

periodic = [1, 1, 0]

trans_field = ffip.Volume_fields_dft(
    center=ffip.Vector3(0, 0, 10*dx),
    size=ffip.Vector3(size.x, size.y, 0), 
    frequency=ft.tolist(), 
    field_component='Ex')

inc_field = ffip.Volume_fields_dft(
    center=ffip.Vector3(0, 0, -10*dx),
    size=ffip.Vector3(size.x, size.y, 0), 
    frequency=ft.tolist(), 
    field_component='Ex')

sources = [ffip.Source(
    function=src_func, 
    center=ffip.Vector3(0, 0, -15*dx), 
    size=ffip.Vector3(size.x, size.y, 0), 
    field_component='Ex') ]

output_filename1 = 'output1.h5'

sim1 = ffip.Simulation(
    size=size, 
    resolution=1/dx,
    pmls=pmls,
    sources=sources,
    fields_output_file=output_filename1,
    fields_output=[inc_field],
    periodic=periodic)

sim1.run(stop_condition=ffip.run_until_time(src_func.end_time * 10))

output_filename2 = 'output2.h5'
sim2 = ffip.Simulation(
    size=size,
    resolution=1/dx,
    pmls=pmls,
    sources=sources,
    geometry=geometry,
    fields_output_file=output_filename2,
    fields_output=[inc_field, trans_field],
    periodic=periodic
)

sim2.run(stop_condition=ffip.run_until_time(src_func.end_time*10))

f1 = h5py.File(output_filename1, 'r')
f2 = h5py.File(output_filename2, 'r')

pts = np.stack(np.meshgrid(ft, 0, 0, 0, indexing='ij'), axis=-1)

inc_interp = inc_field.read(f1, interpolant=1)
ref_interp = inc_field.read(f2, interpolant=1, substract=inc_field.read(f1, interpolant=0))
trn_interp = trans_field.read(f2, interpolant=1)

inc = np.squeeze(inc_interp(pts))
ref = np.squeeze(ref_interp(pts))
trn = np.squeeze(trn_interp(pts))

R = ref / inc
T = trn / inc

plt.figure(1)
plt.plot(l, np.abs(Rt), linewidth=0.5)
plt.plot(l, np.abs(R), '.', markersize=0.5)
plt.title('Reflectance vs Wavelength')

plt.figure(2)
plt.plot(l, np.abs(Tt), linewidth=0.5)
plt.plot(l, np.abs(T), '.', markersize=0.5)
plt.title('Transmittance vs Wavelength')

plt.show()

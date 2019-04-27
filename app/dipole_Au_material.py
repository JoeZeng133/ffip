#%% configure simulation
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

dx = 2e-3/ffip.um_scale
sc = 0.5
dt = sc * dx
l = 500e-3/ffip.um_scale
fs = 1.0/l
ft = 1.0/l

omega = 2 * pi * ft
k0 = omega
eta0 = 1
material = ffip.Au
er = material.get_epsilon(ft)
ur = 1
c = 1 / np.sqrt(er * ur)
eta = np.sqrt(ur / er)
k = omega / c
G = 1

rho = np.linspace(5, 10, 10) * dx
phi = np.linspace(0, 2 * np.pi, 10)
th = np.linspace(0, np.pi, 10)

#indexing = 'ij' equals ndgrid in matlab
Ft, Th, Phi, Rho = np.meshgrid(ft, th, phi, rho, indexing='ij')
Xo, Yo, Zo = ffip.sph2cart(Phi, np.pi / 2 - Th, Rho)

src_func = ffip.Gaussian1(fs, cutoff=4)
t = np.arange(0, src_func.end_time - src_func.start_time, dt)
ref_signal = src_func(t)
ref_signal_fft = np.sum(ref_signal * np.exp(-2j*pi*ft*t), axis=0)

src = ffip.Source(function=src_func, field_component='Ez', amplitude=abs(G))

output_filename = 'output.h5'

sim = ffip.Simulation(
    size=ffip.Vector3(32, 32, 32)*dx, 
    resolution=1/dx,
    courant=sc,
    fields_output_file=output_filename,
    default_material=material,
    pmls=[ffip.PML(thickness=6*dx)],
    sources=[src])

ex_dft = sim.add_dft_fields(size=ffip.Vector3(20, 20, 20) * dx, frequency=[ft], field_component='Ex')
ey_dft = sim.add_dft_fields(size=ffip.Vector3(20, 20, 20) * dx, frequency=[ft], field_component='Ey')
ez_dft = sim.add_dft_fields(size=ffip.Vector3(20, 20, 20) * dx, frequency=[ft], field_component='Ez')

# array broadcasting is done in the reverse direction due to the c-style matrix in python
Erp = 1 * G * dx**3 / (4 * np.pi) * np.exp(-1j * k * Rho) * (2 * eta / Rho**2 + 2 / (1j * omega * er * Rho**3)) * np.cos(Th)
Ethp = 1 * G * dx**3 / (4 * np.pi) * np.exp(-1j * k * Rho) * (1j * omega * ur / Rho + 1 / (1j * omega * er * Rho**3) + 1 * eta / Rho**2) * np.sin(Th)

sim.run(stop_condition=ffip.run_until_time(src_func.end_time * 5))

#%% read fields
fex = ex_dft.get_interpolant(method='nearest',bounds_error=False, fill_value=None)
fey = ey_dft.get_interpolant(method='nearest',bounds_error=False, fill_value=None)
fez = ez_dft.get_interpolant(method='nearest',bounds_error=False, fill_value=None)

pts = np.stack((Ft, Zo, Yo, Xo), axis=-1)
# norm = ref_signal_fft[:, np.newaxis, np.newaxis, np.newaxis]
E = np.stack((fex(pts), fey(pts), fez(pts)), axis=-1) / ref_signal_fft

proj_r = np.stack((np.sin(Th) * np.cos(Phi), np.sin(Th) * np.sin(Phi), np.cos(Th)), axis=-1)
proj_th = np.stack((np.cos(Th) * np.cos(Phi), np.cos(Th) * np.sin(Phi), -np.sin(Th)), axis=-1)

Er = np.sum(E * proj_r, axis=-1)
Eth = np.sum(E * proj_th, axis=-1)


#%% debug
# plt.figure(1)
# plt.plot(np.real(Er.ravel()))
# plt.figure(2)
# plt.plot(np.real(Erp.ravel()))
# plt.show()


#%% plot figures
fig, axes = plt.subplots(1, 2)
axes[0].plot(np.real(Er.ravel()), np.real(Erp.ravel()), '.', markersize=0.5)
axes[0].plot(np.real(Erp.ravel()), np.real(Erp.ravel()), '-', linewidth=0.5)
axes[0].set_title('Re(E_r)')
axes[0].axis('square')

axes[1].plot(np.imag(Er.ravel()), np.imag(Erp.ravel()), '.', markersize=0.5)
axes[1].plot(np.imag(Erp.ravel()), np.imag(Erp.ravel()), '-', linewidth=0.5)
axes[1].set_title('Im(E_r)')
axes[1].axis('square')

fig, axes = plt.subplots(1, 2)
axes[0].plot(np.real(Eth.ravel()), np.real(Ethp.ravel()), '.', markersize=0.5)
axes[0].plot(np.real(Ethp.ravel()), np.real(Ethp.ravel()), '-', linewidth=0.5)
axes[0].set_title('Re(E_th)')
axes[0].axis('square')

axes[1].plot(np.imag(Eth.ravel()), np.imag(Ethp.ravel()), '.', markersize=0.5)
axes[1].plot(np.imag(Ethp.ravel()), np.imag(Ethp.ravel()), '-', linewidth=0.5)
axes[1].set_title('Im(E_th)')
axes[1].axis('square')

plt.show()

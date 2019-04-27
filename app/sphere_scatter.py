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
dx = 2
sc = 0.5
dt = sc * dx

fmin = 1/(1000e-3/ffip.um_scale)
fmax = 1/(200e-3/ffip.um_scale)
fcen = (fmin + fmax) / 2
r=30e-3/ffip.um_scale

material = ffip.Au

size = ffip.Vector3(100, 100, 100)

src_func = ffip.Gaussian1(frequency=fcen)

pmls = [ffip.PML(thickness=10e-3/ffip.um_scale)]

geometry = [ffip.Sphere(radius=r, material=material)]

sources = [ffip.Source(
    function=src_func, 
    center=ffip.Vector3(size.x, 0, -30), 
    size=ffip.Vector3(size.x, size.y, 0), 
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
    center=ffip.Vector3(0, 0, -25),
    size=ffip.Vector3(5, 5, 0), 
    frequency=ft
)

tot_dft = sim2.add_flux_region(
    center=ffip.Vector3(0, 0, -25),
    size=ffip.Vector3(5, 5, 0), 
    frequency=ft
)

trn_dft = sim2.add_flux_region(
    center=ffip.Vector3(0, 0, 25),
    size=ffip.Vector3(5, 5, 0), 
    frequency=ft
)

#%% run simulation
interval = src_func.end_time / 2
sim1.run(stop_condition=ffip.run_until_fields_decay(field_component='Ex', time_interval_examined=interval))
sim2.run(stop_condition=ffip.run_until_fields_decay(field_component='Ex', time_interval_examined=interval))

#%%
omega = 2 * pi * ft
er = material.get_epsilon(ft)
ur = 1
eta = np.sqrt(ur / er)
k = omega * np.sqrt(ur * er)
pts = np.stack(np.meshgrid(ft, 0, 0, 0, indexing='ij'), axis=-1)

inc = inc_dft.get_values(scale=2)
ref = tot_dft.get_values(substract=inc_dft, scale=2)
trn = trn_dft.get_values(scale=2)

R = np.sqrt(-ref / inc)
T = np.sqrt(trn / inc)

#%% plotting

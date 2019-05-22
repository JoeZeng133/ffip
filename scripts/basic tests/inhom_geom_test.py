#%% modules import
import numpy as np
import miepython as mie
import matplotlib.pyplot as plt
import h5py
import ffip
import json
import meep as mp
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#%%
dx = 2
dt = dx * 0.5
lmin = 300
lmax = 1000
nfreq = 50
r = 60
fmin = 1/lmax
fmax = 1/lmin
fcen = (fmin + fmax) / 2
dpml = 20

def density_func(x):
    return 1.0 * (np.sum(x**2, axis=-1) < r**2)

ft = np.linspace(fmin, fmax, nfreq)

material = ffip.Au

size = ffip.Vector3(170, 170, 170)

src_func = ffip.Gaussian1(frequency=fcen)

pmls = [ffip.PML(thickness=dpml)]

vol_size = ffip.Vector3(r*2, r*2, r*2) + 10
geometry = [ffip.Two_Medium_Box(
    size=vol_size, 
    dim=vol_size/dx+1,
    density=density_func,
    medium1=ffip.Au,
    medium2=ffip.Medium(),
    suffix="1" 
    )]

er = material.get_epsilon(frequency=ft)

period = [0, 0, 0]

#%% calculating theorical mie scattering efficiency
m = np.sqrt(er)
x = 2 * pi * r * ft
qext = np.zeros(nfreq)
qabs = np.zeros(nfreq)
qsca = np.zeros(nfreq)

for i in range(nfreq):
    qext[i], qsca[i], qback, g = mie.mie(m[i], x[i])

qabs = qext - qsca
with h5py.File('mie_material.h5', mode='w') as file:
    file.create_dataset('m_real', data=np.real(m))
    file.create_dataset('m_imag', data=np.imag(m))
    file.create_dataset('x', data=x)
    file.create_dataset('l', data=1/ft)
    file.close()
#%%
# plane wave can enter pml areas
sources = [
    ffip.Source(
        function=src_func, 
        center=ffip.Vector3(0, 0, -size.z/2+dpml), 
        size=ffip.Vector3(size.x, size.y, 0), 
        field_component='Ex')]

output_filename1 = 'output1.h5'
output_filename2 = 'output2.h5'

sim1 = ffip.Simulation(
    size=size, 
    resolution=1/dx,
    pmls=pmls,
    sources=sources,
    fields_output_file=output_filename1,
    periodic=period,
    input_file="input1.h5"
)

sim2 = ffip.Simulation(
    size=size,
    resolution=1/dx,
    pmls=pmls,
    sources=sources,
    geometry=geometry,
    fields_output_file=output_filename2,
    periodic=period,
    input_file="input2.h5"
)

# mon_region1 = sim1.add_dft_fields(size=ffip.Vector3(size.x, size.y, size.z), frequency=[fcen], field_component='Ex')
# mon_region2 = sim2.add_dft_fields(size=ffip.Vector3(size.x, size.y, size.z), frequency=[fcen], field_component='Ex')
# interval = src_func.end_time
# sim1.run(stop_condition=ffip.run_until_time(time=interval*2))
# sim2.run(stop_condition=ffip.run_until_time(time=interval*5))
# values = mon_region2.values - mon_region1.values

# with h5py.File('test_output.h5', 'w') as f:
#     f.create_dataset('val_real', data=np.real(values))
#     f.create_dataset('val_imag', data=np.imag(values))
#     f.close()

# input("debug ended")

flux_box_size = ffip.Vector3(r*2, r*2, r*2) + 5

inc_dft = sim1.add_flux_box(
    size=flux_box_size, 
    frequency=ft
)

tot_dft = sim2.add_flux_box(
    size=flux_box_size, 
    frequency=ft
)

#%% run simulation
interval = src_func.end_time / 2
pt1 = ffip.Vector3(0, 0, +flux_box_size.z/2)
pt2 = ffip.Vector3(0, 0, -r)

sim1.run(stop_condition=ffip.run_until_fields_decay(position=pt1, field_component='Ex', time_interval_examined=interval))
sim2.run(stop_condition=ffip.run_until_fields_decay(position=pt2, field_component='Ex', time_interval_examined=interval))

inc_region = inc_dft.get_flux_region(direction='z', side='negative')
inc_power = -inc_region.get_values(scale=1) / (flux_box_size.x * flux_box_size.y) * pi * r**2
absorption = tot_dft.get_values(scale=1) / inc_power * 2
scattering = tot_dft.get_values(scale=1, substract=inc_dft) / inc_power * 2

#%%
plt.subplot(121)
plt.plot(1/ft, -absorption, '.')
plt.plot(1/ft, qabs)
plt.title('absorption')

plt.subplot(122)
plt.plot(1/ft, scattering, '.')
plt.plot(1/ft, qsca)
plt.title('scattering')

plt.show()

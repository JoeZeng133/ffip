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
dpml = 10e-3/ffip.um_scale

l = np.linspace(lmax, lmin, 100)
ft = 1 / l

# total reflection and transmission

size = ffip.Vector3(50e-3/ffip.um_scale, 50e-3/ffip.um_scale, 50e-3/ffip.um_scale)

src_func = ffip.Gaussian1(frequency=fcen)

pmls = [ffip.PML(thickness=dpml)]

def gaussian_beam(sigma, _k, _x0):
    k = np.array([_k.z, _k.y, _k.x])
    x0 = np.array([_x0.z, _x0.y, _x0.x])

    def _gaussian_beam(x):
        return np.exp(1j*2*pi*np.sum(k*(x-x0), axis=-1)-np.sum((x-x0)*(x-x0), axis=-1)/(2*sigma**2))
    return _gaussian_beam

src_size = ffip.Vector3(size.x-dpml*2, size.y-dpml*2) 
dim = src_size/dx + 1

inhom_src = ffip.Inhom_Source(
    function=src_func,
    amplitude=gaussian_beam(2, ffip.Vector3(z=1), ffip.Vector3(z=-size.z/2+dpml)),
    dim=dim,
    center=ffip.Vector3(z=-size.z/2+dpml),
    size=src_size,
    frequency=fcen,
    field_component='Ex',
    suffix="1"
)

amp = np.squeeze(inhom_src.amplitude)
plt.imshow(np.abs(amp))
plt.show()

sources = [inhom_src]

periodic = [0, 0, 0]

sim1 = ffip.Simulation(
    size=size, 
    resolution=1/dx,
    pmls=pmls,
    sources=sources,
    fields_output_file='gbeam_output.h5',
    input_file='gbeam_input.h5',
    periodic=periodic
)

xyplane = sim1.add_dft_fields(size=ffip.Vector3(size.x, size.y, 0), frequency=[fcen], field_component='Ex')
xzplane = sim1.add_dft_fields(size=ffip.Vector3(size.x, 0, size.z), frequency=[fcen], field_component='Ex')

#%%
sim1.run(stop_condition=ffip.run_until_time(time=src_func.end_time*4))

#%%

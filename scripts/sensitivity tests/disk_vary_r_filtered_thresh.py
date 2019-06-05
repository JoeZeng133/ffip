#%% modules import
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interpn, RegularGridInterpolator
from scipy.signal import convolve
import h5py
import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import linprog

def disk(r0=30, a=2):
    def _disk(x):
        r = np.sqrt(np.sum(x**2, axis=-1))
        return -1/np.pi * np.arctan(a * (r - r0)) + 0.5
    
    return _disk

def disk_der(r0=30, a=2):
    def _disk(x):
        r = np.sqrt(np.sum(x**2, axis=-1))
        return a / np.pi * 1 / (1 + (a * (r - r0))**2)
    
    return _disk

#%% original simulation
prefix = 'disk_vary_r_filtered_thresh_'
dx = 4
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(40, 40, 40) * dx + dpml * 2
fsrc = 1 / 500
fcen = 1 / 800

m1 = ffip.Au
m2 = ffip.Medium()

e1 = m1.get_epsilon(fcen)
e2 = m2.get_epsilon(fcen)

disk_h = 80

src_func = ffip.Gaussian1(fsrc, start_time=0.5/fcen)

ref_t = np.arange(src_func.start_time, src_func.end_time, dt)
ref_fft = np.sum(np.exp(-2j*pi*fcen*ref_t) * src_func(ref_t))

pmls = [ffip.PML(
    thickness=dpml
)]

sources = [ffip.Source(
    function=src_func,
    center=ffip.Vector3(z=-sim_size.z/2+dpml),
    size=ffip.Vector3(sim_size.x, sim_size.y, 0),
    field_component='Ex')
]

adj_source_size = ffip.Vector3(x=120, y=120)
adj_source_center = ffip.Vector3(z=50)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))

geom_size = ffip.Vector3(120, 120, disk_h)
geom_dim = (geom_size/dx+1).round()
geom_center = ffip.Vector3()

stop_condition = ffip.run_until_dft(
    center=geom_center,
    size=geom_size,
    field_component='Ex',
    time_interval_examined=4/fsrc,
    var=5e-3,
    frequency=fcen
)

e_vals = np.zeros(adj_source_shape)

#%% optimization
sim_forward = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    sources=sources,
    pmls=pmls,
    input_file=prefix+'forward_input.h5',
    fields_output_file=prefix+'forward_output.h5',
    progress_interval=1e10
)

sim_adjoint = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    input_file=prefix+'adjoint_input.h5',
    fields_output_file=prefix+'adjoint_output.h5',
    progress_interval=1e10
)

adj_src = ffip.Adjoint_Source(
    adjoint_simulation=sim_adjoint,
    forward_simulation=sim_forward,
    function=src_func,
    frequency=fcen,
    center=adj_source_center,
    size=adj_source_size,
    dim=adj_source_dim,
    functionals=[
        ['|E|', e_vals]
    ]
)

adj_vol = ffip.Adjoint_Volume(
    adjoint_simulation=sim_adjoint,
    forward_simulation=sim_forward,
    frequency=fcen,
    center=geom_center,
    size=geom_size,
    dim=geom_dim,
    density=None,
    medium1=m1,
    medium2=m2,
    norm=ref_fft
)

def get_gaussian_filter(filter_r = 5):

    filter_sigma = filter_r / 2
    filter_x = np.arange(-filter_r, filter_r+1)
    filter_X, filter_Y = np.meshgrid(filter_x, filter_x, indexing='ij')
    filter = np.exp(-0.5 * (filter_X**2 + filter_Y**2) / filter_sigma)
    filter = filter / np.sum(filter.ravel())

    return filter

def filt(x, filter):
    return convolve(x, filter, mode='same')

def filt_transpose(xp, filter):
    tp = np.flip(filter, axis=None)
    return convolve(xp, tp, mode='same')

def thresh(x, beta, eta):
    return (np.tanh(beta * eta) + np.tanh(beta * (x - eta))) / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))

def thresh_transpose(xp, beta, eta):
    return beta * (1 / np.cosh(beta * (xp - eta)))**2 / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))

filter = get_gaussian_filter(3)
rho_shape = adj_vol.density.shape
geom_pts = np.stack(np.meshgrid(adj_vol.y, adj_vol.x, indexing='ij'), axis=-1)
r_list = np.linspace(40, 50, 10)
f_list = []
se1s = []
se2s = []
beta = 5
eta = 0.5

for r in r_list:

    rho_der = disk_der(r0=r)(geom_pts)
    rho = disk(r0=r)(geom_pts)

    # applying filter and threshold projection
    rho = filt(rho, filter)
    rho = thresh(rho, beta=beta, eta=eta)

    plt.figure(1)
    plt.imshow(rho)
    plt.title('r=%f' % r)
    plt.pause(1)

    adj_vol.density = np.ones(rho_shape) * rho

    sim_forward.run(
        # skip=True, pop=True,
        stop_condition=stop_condition, 
        np=3
    )

    f1 =  adj_src.eval_functionals_and_set_sources()
    f_list.append(f1)

    plt.figure(2)
    plt.plot(r, f1, '.')
    plt.pause(1)

    sim_adjoint.run(
        # skip=True, pop=True,
        stop_condition=stop_condition, 
        np=3
    )

    se1 = np.sum(adj_vol.get_sensitivity(), axis=0)
    se2 = np.sum(adj_vol.get_sensitivity2(), axis=0)

    # applying transpose
    se1 = thresh_transpose(se1, beta=beta, eta=eta)
    se1 = filt_transpose(se1, filter)

    se2 = thresh_transpose(se2, beta=beta, eta=eta)
    se2 = filt_transpose(se2, filter)

    se1s.append(np.sum((se1 * rho_der).ravel()))
    se2s.append(np.sum((se2 * rho_der).ravel()))


with h5py.File('tmp2.h5', 'w') as file:
    file.create_dataset('r', data=r_list)
    file.create_dataset('se1', data=np.array(se1s))
    file.create_dataset('se2', data=np.array(se2s))
    file.create_dataset('f', data=np.array(f_list))

plt.show()
input('end')
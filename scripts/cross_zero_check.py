#%% modules import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def disk(r=30):
    def _disk(x):
        return 1.0 * (np.sum(x**2, axis=-1) < r**2)
    
    return _disk

def arb(x):
    return 0.1 + 0.8 * np.abs((np.sin(x[..., 1]/dx/1.5) * np.sin(x[..., 2]/dx/1.5)))

def arb1(x):
    return 0.05 * np.ones(x.shape[:-1])

#%% original simulation
dx = 2
dt = 0.5 * dx
sim_size = ffip.Vector3(50, 50, 50) * dx
dpml = 8 * dx
fsrc = 1 / 500
fcen = 1 / 800

m1 = ffip.Au
m2 = ffip.Medium()

e1 = m1.get_epsilon(fcen)
e2 = m2.get_epsilon(fcen)

rho_cross = -np.real(e2) / np.real(e1 - e2)

print('rho_cross=', rho_cross)

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

adj_source_size = ffip.Vector3(x=40, y=40)
adj_source_center = ffip.Vector3(z=20)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))

geom_size = ffip.Vector3(30, 30, 30)
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

ex_vals = np.zeros(adj_source_shape)
ey_vals = np.zeros(adj_source_shape)
ez_vals = np.zeros(adj_source_shape)

#%% optimization
sim_forward = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    sources=sources,
    pmls=pmls,
    input_file='forward_input.h5',
    fields_output_file='forwad_output.h5',
    progress_interval=1e10
)

sim_adjoint = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    input_file='adjoint_input.h5',
    fields_output_file='adjoint_output.h5',
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
        ['|E|', ex_vals]
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

rho_shape = adj_vol.density.shape
rho_list = np.linspace(0, 0.2, 100)
f_list = []
se1s = []
se2s = []

for rho in rho_list:
    adj_vol.density = np.ones(rho_shape) * rho

    sim_forward.run(stop_condition=stop_condition, np=16)

    f1 =  adj_src.eval_functionals_and_set_sources()
    f_list.append(f1)

    plt.plot(rho, f1, '.')
    plt.pause(1)

    sim_adjoint.run(stop_condition=stop_condition, np=16)

    se1s.append(np.sum(adj_vol.get_sensitivity().ravel()))
    se2s.append(np.sum(adj_vol.get_sensitivity2().ravel()))


with h5py.File('tmp2.h5', 'w') as file:
    file.create_dataset('rho', data=rho_list)
    file.create_dataset('se1', data=np.array(se1s))
    file.create_dataset('se2', data=np.array(se2s))
    file.create_dataset('f', data=np.array(f_list))

plt.show()
input('end')

#%% modules import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import subprocess

#%% original simulation
prefix='block_vary_perm'
num_cores = 35
dx = 2
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(100, 100, 100) + dpml * 2
fsrc = 1 / 500
fcen = 1 / 800

m1 = ffip.Au
m2 = ffip.Medium()

e1 = m1.get_epsilon(fcen)
e2 = m2.get_epsilon(fcen)

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

geom_size = ffip.Vector3(30, 30, 30)
geom_dim = (geom_size/dx+1).round()
geom_center = ffip.Vector3()

stop_condition = ffip.run_until_dft(
    center=geom_center,
    size=geom_size,
    field_component='Ex',
    time_interval_examined=4/fsrc,
    var=1e-3,
    frequency=fcen
)


#%% forward and adjoint set up
sim_forward = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    sources=sources,
    pmls=pmls,
    input_file=prefix + 'forward_input.h5',
    fields_output_file=prefix + 'forwad_output.h5',
    progress_interval=20
)

sim_adjoint = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    input_file=prefix + 'adjoint_input.h5',
    fields_output_file=prefix + 'adjoint_output.h5',
    progress_interval=20
)

box_flux = sim_forward.add_flux_box(
    center=ffip.Vector3(),
    size=geom_size + 10,
    frequency=[fcen]
)

adj_src = ffip.Adjoint_Flux(
    adjoint_simulation=sim_adjoint,
    forward_simulation=sim_forward,
    function=src_func,
    frequency=fcen,
    fluxes=box_flux.flux_regions
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
rho_list = np.linspace(0, 0.25, 1001)
f_list = []
se1s = []
se2s = []

for rho in rho_list:
    adj_vol.density = np.ones(rho_shape) * rho

    print('running forward simulation')
    sim_forward.run(
        stdout=subprocess.DEVNULL,
        stop_condition=stop_condition, 
        np=num_cores
    )

    f =  adj_src.eval_functionals_and_set_sources()
    f_list.append(f)

    print('rho=%f -> f=%f' % (rho, f))

    print('running adjoint simulation')

    sim_adjoint.run(
        stdout=subprocess.DEVNULL,
        stop_condition=stop_condition,
        np=num_cores
    )

    se1s.append(np.sum(adj_vol.get_sensitivity().ravel()))
    se2s.append(np.sum(adj_vol.get_sensitivity2().ravel()))

    with h5py.File(prefix + 'result.h5', 'w') as file:
        file.create_dataset('rho', data=rho_list)
        file.create_dataset('se1', data=np.array(se1s))
        file.create_dataset('se2', data=np.array(se2s))
        file.create_dataset('f', data=np.array(f_list))




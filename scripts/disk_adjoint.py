#%% modules import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ffip
import json
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import linprog

def disk(r=30):
    def _disk(x):
        return 1.0 * ((x[..., 1]**2 + x[..., 2]**2) < r**2)
    
    return _disk


#%% original simulation
dx = 2
dt = 0.5 * dx
sim_size = ffip.Vector3(66, 66, 66) * dx
dpml = 8 * dx
fsrc = 1 / 800
fcen = 1 / 800
disk_r = 30
disk_h = 40

m1 = ffip.Au
m2 = ffip.Medium()

e1 = m1.get_epsilon(fcen)
e2 = m2.get_epsilon(fcen)

e1d = m1.get_dis_epsilon(fcen, dt)
e2d = m2.get_dis_epsilon(fcen, dt)

print('continuous e1=', e1, ',e2=', e2)
print('digital e1=', e1d, ',e2=', e2d)

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

adj_source_size = ffip.Vector3(x=80, y=80)
adj_source_center = ffip.Vector3(z=30)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))

geom_size = ffip.Vector3(disk_r * 2 + 10, disk_r * 2 + 10, disk_h)
geom_dim = (geom_size/dx/2+1).round()
geom_center = ffip.Vector3()

# stop_condition = ffip.run_until_fields_decay(
#     position=ffip.Vector3(z=15),
#     field_component='Ex',
#     time_interval_examined=1/fsrc,
#     decayed_by=5e-4
# )

# stop_condition = ffip.run_until_time(
#     time=dt*20000
# )

stop_condition = ffip.run_until_dft(
    center=geom_center,
    size=geom_size,
    field_component='Ex',
    time_interval_examined=1/fsrc,
    var=1e-2,
    frequency=fcen
)

geom0 = [ffip.Two_Medium_Box(
    size=geom_size,
    center=geom_center,
    dim=geom_dim,
    density=disk(r=disk_r),
    medium1=m1,
    medium2=m2,
    suffix='0'
)]

sim0 = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    geometry=geom0,
    sources=sources,
    progress_interval=10,
    input_file='sim0_input.h5',
    fields_output_file='sim0_output.h5'
)

ex_dft = sim0.add_dft_fields(
    center=adj_source_center,
    size=adj_source_size,
    frequency=[fcen],
    field_component='Ex'
)

ey_dft = sim0.add_dft_fields(
    center=adj_source_center,
    size=adj_source_size,
    frequency=[fcen],
    field_component='Ey',
)

ez_dft = sim0.add_dft_fields(
    center=adj_source_center,
    size=adj_source_size,
    frequency=[fcen],
    field_component='Ez'
)

# add for stop condition
sim0.add_dft_fields(
    center=geom_center,
    size=geom_size,
    frequency=[fcen],
    field_component='Ex'   
)

sim0.run(stop_condition=stop_condition, np=1, skip=False)

adj_x, adj_y, adj_z = ffip.getgrid(adj_source_center, adj_source_size, adj_source_dim)

ex_vals = ex_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
ey_vals = ey_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
ez_vals = ez_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)

plt.figure(1)
plt.imshow(geom0[0].density[0, ...], vmin=0, vmax=1)
plt.colorbar()
plt.title('Original Shape')
# plt.draw()
plt.pause(1)

# %% optimization

sim_forward = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    sources=sources,
    pmls=pmls,
    input_file='forward_input.h5',
    fields_output_file='forwad_output.h5',
    progress_interval=20
)

sim_adjoint = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    input_file='adjoint_input.h5',
    fields_output_file='adjoint_output.h5',
    progress_interval=20
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
        ['Ex', ex_vals],
        ['Ey', ey_vals],
        ['Ez', ez_vals]
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

max_itr = 10
itr = 0
funcs = []
funcs_exp = []
rhos = []
step = 2e-2
decay = 0.95
decay_after_step = 10

rho = np.array(geom0[0].density[0,...]) * 0.8
# rho = rho + np.random.random(rho.shape) * 0.1
# rho[rho < 0] = 0
# rho[rho > 1] = 1

# with h5py.File('disk_adjoint_rho1.h5', 'r') as file:
#     rho = file['rho 0'][:]
#     file.close()

file = h5py.File('disk_adjoint_rho.h5', 'w')

while itr < max_itr:
    print('##################################################')
    print('Iteration ', itr)

    plt.figure(itr+2)
    plt.imshow(rho, vmin=0, vmax=1)
    plt.colorbar()
    plt.title('rho at iteration %d' % itr)
    # plt.draw()
    plt.pause(1)

    file.create_dataset('rho %d' % itr, data=rho)

    # broadcast to 3d rho
    adj_vol.density = np.ones(adj_vol.shape) * rho

    sim_forward.run(stop_condition=stop_condition, np=12)

    funcs.append(adj_src.eval_functionals_and_set_sources())

    print('Objective function=', funcs[itr])
    if itr > 0:
        print('Diff=', funcs[itr] - funcs[itr-1])

    sim_adjoint.run(stop_condition=stop_condition, np=12)

    se = adj_vol.get_sensitivity()

    se_2d = np.sum(se, axis=0)

    # optimization, linear programming
    lb = np.maximum(-rho.ravel(), -step)
    ub = np.minimum(1 - rho.ravel(), step)
    bounds = np.stack((lb, ub), axis=-1)

    res = linprog(c=se_2d.ravel(), bounds=bounds)

    diff = np.sum(res.x * se_2d.ravel())

    print('Expected diff=', diff)

    rho = rho + np.reshape(res.x, rho.shape)

    

    itr += 1
    if itr > decay_after_step:
        step *= decay
# %%
file.close()
plt.show()
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

#%% original simulation
dx = 2
dt = 0.5 * dx
sim_size = ffip.Vector3(50, 50, 50) * dx
dpml = 8 * dx
fsrc = 1 / 500
fcen = 1 / 800

m1 = ffip.Au
m2 = ffip.Medium()

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

# stop_condition = ffip.run_until_fields_decay(
#     position=ffip.Vector3(z=15),
#     field_component='Ex',
#     time_interval_examined=1/fsrc,
#     decayed_by=5e-4
# )

stop_condition = ffip.run_until_time(
    time=dt*10000
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
    center=ffip.Vector3(),
    size=geom_size,
    dim=geom_dim,
    density=None,
    medium1=m1,
    medium2=m2,
    norm=ref_fft
)

# adj_vol.density = adj_vol.density + 0.5
plt.imshow(adj_vol.density[0, :, :])
plt.show()

sim_forward.run(stop_condition=stop_condition, np=10)

f1 =  adj_src.eval_functionals_and_set_sources()
print('Objective function is evaluated at', f1)

sim_adjoint.run(stop_condition=stop_condition, np=10)

se1 = adj_vol.get_sensitivity()
se2 = adj_vol.get_sensitivity2()

pert = np.random.random(adj_vol.density.shape) * 2e-2
adj_vol.density = adj_vol.density + pert

diff1 = np.sum((se1 * pert).ravel())
diff2 = np.sum((se2 * pert).ravel())

print('exp diff1=', diff1)
print('exp diff2=', diff2)

sim_forward.input_file = 'change_input.h5'
sim_forward.fields_output_file = 'change_output.h5'

sim_forward.run(stop_condition=stop_condition, np=10)

f2 =  adj_src.eval_functionals_and_set_sources()

print('Changed Objective function is evaluated at', f2)
print('actual diff=', f2 - f1)
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

#%% original simulation
dx = 2
dt = 0.5 * dx
fcen = 1 / 500
sim_size = ffip.Vector3(130, 130, 100)
dpml = 15
fsrc = 1 / 500
fcen = 1 / 500
disk_r = 30
disk_h = 40

m1 = ffip.Au
m2 = ffip.Medium()

disk_size = ffip.Vector3(disk_r*2+10, disk_r*2+10, disk_h)
src_func = ffip.Gaussian1(fsrc, start_time=0.5/fcen)

pmls = [ffip.PML(
    thickness=dpml
)]

geom0 = ffip.Two_Medium_Box(
    size=disk_size, 
    center=ffip.Vector3(),
    dim=disk_size/dx+1, 
    density=disk(disk_r),
    medium1=m1,
    medium2=m2,
    suffix='0'
)

geometry = [geom0]

sources = [ffip.Source(
    function=src_func,
    center=ffip.Vector3(z=-sim_size.z/2+dpml),
    size=ffip.Vector3(sim_size.x, sim_size.y, 0),
    field_component='Ex')
]

sim0 = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    geometry=geometry,
    sources=sources,
    input_file='sim0_input.h5',
    fields_output_file='sim0_out.h5'
)

adj_source_size = ffip.Vector3(disk_r*2, disk_r*2)
adj_source_center = ffip.Vector3(z=disk_h/2+10)

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
    field_component='Ey'
)

ez_dft = sim0.add_dft_fields(
    center=adj_source_center, 
    size=adj_source_size,
    frequency=[fcen],
    field_component='Ez'
)

stop_condition = ffip.run_until_fields_decay(
    position=ffip.Vector3(z=disk_h/2+10), 
    field_component='Ex', 
    time_interval_examined=1/fsrc,
    decayed_by=1e-3)

sim0.run(stop_condition=stop_condition)

#%% retrieve objective values
adj_source_dim = (adj_source_size/dx).round() + 1
adj_source_p1 = adj_source_center - adj_source_size/2
adj_source_x = np.linspace(adj_source_p1.x, adj_source_p1.x + adj_source_size.x, adj_source_dim.x)
adj_source_y = np.linspace(adj_source_p1.y, adj_source_p1.y + adj_source_size.y, adj_source_dim.y)
adj_source_z = np.linspace(adj_source_p1.z, adj_source_p1.z + adj_source_size.z, adj_source_dim.z)

ex_vals = ex_dft(fcen, adj_source_z, adj_source_y, adj_source_x)
ey_vals = ey_dft(fcen, adj_source_z, adj_source_y, adj_source_x)
ez_vals = ez_dft(fcen, adj_source_z, adj_source_y, adj_source_x)

#%% optimization
sim_forward = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    sources=sources,
    pmls=pmls,
    input_file='forward_input.h5',
    fields_output_file='forwad_output.h5'
)

sim_adjoint = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    input_file='adjoint_input.h5',
    fields_output_file='adjoint_output.h5'
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
    center=ffip.Vector3(),
    size=disk_size,
    dim=geom0.dimension,
    density=None,
    medium1=m1,
    medium2=m2
)

density = adj_vol.density

sim_forward.run(stop_condition=stop_condition)

f1 =  adj_src.eval_functionals_and_set_sources()
print('Objective function is evaluated at', f1)

sim_adjoint.run(stop_condition=stop_condition)

sensitivity = adj_vol.get_sensitivity()

adj_vol.density = adj_vol.density + 5e-2

diff = np.sum(sensitivity.ravel() * 5e-2)

sim_forward.run(stop_condition=stop_condition)

f2 =  adj_src.eval_functionals_and_set_sources()

print('Changed Objective function is evaluated at', f2)
print('exp diff=', diff, ', actual diff=', f2 - f1)

















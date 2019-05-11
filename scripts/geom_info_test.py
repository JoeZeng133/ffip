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
r = 30
dpml = 20

material = ffip.Au

size = ffip.Vector3(120, 120, 120)

pmls = [ffip.PML(thickness=dpml)]

# geometry = [ffip.Sphere(radius=r, material=material)]
geometry = [ffip.Box(size=ffip.Vector3(r, r, r))]

sim1 = ffip.Simulation(
    size=size, 
    resolution=1/dx,
    pmls=pmls,
    geometry=geometry,
    fields_output_file='geom_info.h5'
)

ex_pts = sim1.request_geom_info(geometry[0], field_component='Ex')
ey_pts = sim1.request_geom_info(geometry[0], field_component='Ey')
ez_pts = sim1.request_geom_info(geometry[0], field_component='Ez')

sim1.run(stop_condition=ffip.run_until_time(time=0))



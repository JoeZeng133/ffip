import meep as mp
from meep.materials import Au
import ffip
from math import pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

def meep_test():
    dx = 2
    sc = 0.5
    dt = sc * dx
    fmax = 1/400
    fmin = 1/800
    nf = 20
    resolution = 1/dx

    # parameters
    cell_size = mp.Vector3(34, 34, 34) * dx
    sources = [mp.Source(mp.GaussianSource(frequency=(fmin + fmax) / 2, fwidth=fmax - fmin), component=mp.Ez, center=mp.Vector3())]
    pml_layers = [mp.PML(6*dx)]
    geometry = [mp.Block(size=mp.Vector3(1e9, 1e9, 1e9), center=mp.Vector3(), material=Au)]
    sim = mp.Simulation(cell_size=cell_size,
        dimensions=3,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution)

    dft_vol = mp.Volume(center=mp.Vector3(), size=cell_size, dims=3)
    dft_fields = sim.add_dft_fields(components=[mp.Ex, mp.Ey, mp.Ez], 
        freq_min=fmin, 
        freq_max=fmax, 
        nfreq=nf,
        where=dft_vol)
    
    sim.run(until=5000*dt)
    field = []
    for i in range(20):
        field.append(sim.get_dft_array(dft_fields, component=mp.Ex, num_freq=0))
    return field

def ffip_test():
    dx = 2
    sc = 0.5
    dt = sc * dx
    fmax = 1/400
    fmin = 1/800
    fs = (fmax + fmin) / 2
    nf = 20
    ft = np.linspace(fmin, fmax, nf)
    resolution = 1/dx

    src_func = ffip.Gaussian1(fs, cutoff=4)
    cell_size = ffip.Vector3(34, 34, 34) * dx
    sources = [ffip.Source(function=src_func, field_component='Ez', amplitude=1)]
    pml_layers = [ffip.PML(thickness=6*dx)]
    default_material = ffip.Au
    ex_dft = ffip.Volume_fields_dft(size=cell_size, frequency=ft.tolist(), field_component='Ex')

    output_filename = 'output1.h5'
    sim = ffip.Simulation(size=cell_size,
        resolution=resolution,
        pmls=pml_layers,
        default_material=default_material,
        fields_output=ex_dft,
        sources=sources,
        fields_output_file=output_filename
        )
    
    sim.run(stop_condition=ffip.run_until_time(5000*dt))
    input('hit anything when finished')
    field = ex_dft.read(fo, interpolant=0)
    return field




import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from meep.materials import Au
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import ffip

def meep_sim():
    dx = 2e-3
    resolution = 1/dx
    lmin = 200e-3
    lmax = 1000e-3
    fmin = 1/lmax
    fmax = 1/lmin
    nfreq = 100

    fcen = (fmin + fmax) / 2
    df = fmax - fmin
    d = 40e-3

    cell_size = mp.Vector3(10e-3, 10e-3, 100e-3)
    pml_layers = [mp.PML(20e-3, direction=mp.Z, side=mp.High), mp.PML(20e-3, direction=mp.Z, side=mp.Low)]
    src = [mp.Source(
        mp.GaussianSource(fcen, fwidth=df),
        component=mp.Ex,
        center=mp.Vector3(0, 0, -30e-3),
        size=mp.Vector3(cell_size.x, cell_size.y, 0))]
    geometry = [mp.Block(mp.Vector3(mp.inf, mp.inf, d), material=Au, center=mp.Vector3())]

    sim = mp.Simulation(
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=src,
        k_point=mp.Vector3(),
        resolution=resolution)

    before = mp.FluxRegion(center=mp.Vector3(0, 0, -25e-3), size=mp.Vector3(10e-3, 10e-3, 0))
    after  = mp.FluxRegion(center=mp.Vector3(0, 0, 25e-3), size=mp.Vector3(10e-3, 10e-3, 0))

    inc = sim.add_flux(fcen, df, nfreq, before)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), decay_by=1e-4))

    inc_flux = mp.get_fluxes(inc)
    inc_flux_data = sim.get_flux_data(inc)
    sim.reset_meep()

    sim = mp.Simulation(
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=src,
        geometry=geometry,
        k_point=mp.Vector3(),
        resolution=resolution
    )

    trn = sim.add_flux(fcen, df, nfreq, after)
    ref = sim.add_flux(fcen, df, nfreq, before)

    sim.load_minus_flux_data(ref, inc_flux_data)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), decay_by=1e-4))

    trn_flux = mp.get_fluxes(trn)
    ref_flux = mp.get_fluxes(ref)
    freqs = np.asarray(mp.get_flux_freqs(trn))
    l = 1/freqs

    Rm = np.sqrt(-np.asarray(ref_flux) / np.asarray(inc_flux))
    Tm = np.sqrt(np.asarray(trn_flux) / np.asarray(inc_flux))

    return Rm, Tm

d = 40
f = freqs / ffip.um_scale
omega = 2 * pi * f
material = ffip.Au
er = material.get_epsilon(f)
ur = 1
eta = np.sqrt(ur / er)
k = omega * np.sqrt(ur * er)

R1 = (eta - 1) / (eta + 1)
T1 = R1 + 1
R2 = (1 - eta) / (eta + 1)
T2 = R2 + 1

# total reflection and transmission
Rt = (R1 + R2 * np.exp(-2j * k * d)) / (1 + R1 * R2 * np.exp(-2j * k * d))
Tt = (T1 * T2 * np.exp(-1j * k * d)) / (1 + R1 * R2 * np.exp(-2j * k * d))



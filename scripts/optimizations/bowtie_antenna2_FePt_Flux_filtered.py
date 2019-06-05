import ffip
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import minimize, Bounds
from scipy.signal import convolve
import subprocess

#%% setting up target geometry
def get_bowtie_den():
    a = 20
    b = 40
    c = 200
    d = 200

    tri1 = ffip.Polygon([(a/2, b/2), (c/2, d/2), (c/2, -d/2), (a/2, -b/2)])
    tri2 = ffip.Polygon([(-a/2, b/2), (-c/2, d/2), (-c/2, -d/2), (-a/2, -b/2)])

    res = tri1.union(tri2)
    density = ffip.planar_polygon(res)

    return density

density_fun = get_bowtie_den()

#%% setting up basic configurations
prefix = 'bowtie_antenna2_FePt_Flux_filtered_'
dx = 4
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(300, 300, 120) + dpml * 2
fsrc = 1 / 500
fcen = 1 / 800

m1 = ffip.Au
m2 = ffip.Medium()
fept = ffip.FePt(fcen)

src_func = ffip.Gaussian1(fsrc, start_time=0.5/fcen)
ref_t = np.arange(src_func.start_time, src_func.end_time, dt)
ref_fft = np.sum(np.exp(-2j*np.pi*fcen*ref_t) * src_func(ref_t))

pmls = [ffip.PML(
    thickness=dpml
)]

sources = [ffip.Source(
    function=src_func,
    center=ffip.Vector3(z=-sim_size.z/2+dpml),
    size=ffip.Vector3(sim_size.x, sim_size.y, 0),
    field_component='Ex')
]

adj_source_size = ffip.Vector3(x=50, y=50, z=10)
adj_source_center = ffip.Vector3(z=35)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))

geom_size = ffip.Vector3(200, 200, 60)
geom_dim = (geom_size/dx+1).round()
geom_center = ffip.Vector3(z=-10)

bowtie = ffip.Two_Medium_Box(
    size=geom_size,
    center=geom_center,
    dim=geom_dim,
    density=density_fun,
    medium1=m1,
    medium2=m2
)

bowtie_x = bowtie.x
bowtie_y = bowtie.y
bowtie_dx = bowtie_x[1] - bowtie_x[0]
bowtie_dy = bowtie_y[1] - bowtie_y[0]
extent = (bowtie_x[0] - bowtie_dx/2, bowtie_x[-1] - bowtie_dx/2, bowtie_y[0] - bowtie_dy/2, bowtie_y[-1] - bowtie_dy/2)
bowtie_density = bowtie.density[0, ...]
plt.imshow(bowtie_density, origin='lower', extent=extent)
plt.show()

medium_layer = ffip.Box(
    size=ffip.Vector3(sim_size.x, sim_size.y, 10),
    center=ffip.Vector3(z=35),
    material=fept
)

geometry = [bowtie]

stop_condition = ffip.run_until_dft(
    center=geom_center,
    size=geom_size,
    field_component='Ex',
    time_interval_examined=4/fsrc,
    var=5e-3,
    frequency=fcen
)

#%% adjoint setups
sim_forward = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    sources=sources,
    pmls=pmls,
    input_file=prefix + 'forward_input.h5',
    fields_output_file=prefix + 'forward_output.h5',
    progress_interval=20
)

sim_adjoint = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    input_file=prefix+'adjoint_input.h5',
    fields_output_file=prefix+'adjoint_output.h5',
    progress_interval=20
)

box_flux = sim_forward.add_flux_box(
    center=adj_source_center,
    size=adj_source_size,
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

def get_gaussian_filter(filter_r = 5):

    filter_sigma = filter_r / 2
    filter_x = np.arange(-filter_r, filter_r+1)
    filter_X, filter_Y = np.meshgrid(filter_x, filter_x, indexing='ij')
    filter = np.exp(-0.5 * (filter_X**2 + filter_Y**2) / filter_sigma)
    filter = filter / np.sum(filter.ravel())

    return filter

filter = get_gaussian_filter(5)
filter_transpose = np.flip(filter, axis=None)

scale = 1
rho_list = []
f_list = []

def fun(rho):
    tmp = np.reshape(rho / scale, adj_vol.shape[1:])
    tmp = convolve(tmp, filter, mode='same')

    adj_vol.density = np.ones(adj_vol.shape) * tmp
    print('running forward calculation')
    sim_forward.run(
        stop_condition=stop_condition, 
        np=3,
        stdout=subprocess.DEVNULL
        # skip=True,
        # pop=True
    )
    
    return adj_src.eval_functionals_and_set_sources()

def fprime(rho):
    print('running adjoint calculation')
    sim_adjoint.run(
        stop_condition=stop_condition,
        np=3,
        stdout=subprocess.DEVNULL
        # skip=True,
        # pop=True
    )

    se = np.sum(adj_vol.get_sensitivity() / scale, axis=0)
    se = convolve(se, filter_transpose, mode='same')

    return se.ravel()

def show(rho):
    rho_t = np.reshape(rho / scale, adj_vol.shape[1:])

    rho_list.append(np.array(rho_t))

    f_list.append(adj_src.eval_functionals_and_set_sources())

    plt.figure(1)
    plt.imshow(rho_t, origin='lower', vmin=0, vmax=1)
    plt.pause(3)

#%% sensitivity test
def sensitivity_test():
    rho0 = np.zeros(adj_vol.shape[1:]).flatten() + 0.5
    f1 = fun(rho0)
    print('f1=', f1)
    print('exp(f2 - f1)=', np.sum(fprime(rho0) * 0.01))
    sim_forward.fields_output_file = prefix+'forward_pert_output.h5'
    sim_forward.input_file = prefix+'forward_pert_input.h5'
    f2 = fun(rho0 + 0.01)
    print('f2=', f2)
    print('f2-f1=', f2 - f1)
    print("end")

sensitivity_test()

#%%
def optimize():
    rho0 = np.zeros(adj_vol.shape[1:]).flatten()
    res = minimize(
        fun,
        rho0,
        method='L-BFGS-B',
        jac=fprime,
        bounds=Bounds(0, scale),
        callback=show,
        options={'maxiter' : 20, 'disp' : True}
    )

    print(res)

    with h5py.File(prefix + 'result.h5', 'w') as file:
        # save e rho
        file.create_dataset('rho target', data=bowtie_density)
        file.attrs['rho target x'] = bowtie_x
        file.attrs['rho target y'] = bowtie_y

        # save funs
        file.create_dataset('fun', data=np.array(f_list))

        # save rhos and es
        for i in range(len(rho_list)):
             file.create_dataset('rho %d' % i, data=rho_list[i])

# optimize()
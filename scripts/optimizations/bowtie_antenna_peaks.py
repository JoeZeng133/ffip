import ffip
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import minimize, Bounds
import subprocess

#%% setting up target geometry
def get_bowtie_den():
    a = 30
    b = 60
    c = 400
    d = 400

    tri1 = ffip.Polygon([(a/2, b/2), (c/2, d/2), (c/2, -d/2), (a/2, -b/2)])
    tri2 = ffip.Polygon([(-a/2, b/2), (-c/2, d/2), (-c/2, -d/2), (-a/2, -b/2)])

    res = tri1.union(tri2)
    density = ffip.planar_polygon(res)

    return density

density_fun = get_bowtie_den()

#%% setting up basic configurations
prefix = 'bowtie_antenna_peaks_'
dx = 5
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(500, 500, 150) + dpml * 2
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

adj_source_size = ffip.Vector3(x=200, y=200)
adj_source_center = ffip.Vector3(z=50)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))

geom_size = ffip.Vector3(400, 400, 90)
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
    center=ffip.Vector3(z=50),
    material=fept
)

geometry = [bowtie, medium_layer]

stop_condition = ffip.run_until_dft(
    center=geom_center,
    size=geom_size,
    field_component='Ex',
    time_interval_examined=4/fsrc,
    var=5e-3,
    frequency=fcen
)

#%% run and get target fields
def get_target_fields():
    sim0 = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        pmls=pmls,
        geometry=geometry,
        sources=sources,
        progress_interval=20,
        input_file=prefix+'input.h5',
        fields_output_file=prefix+'output.h5'
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
        field_component='Ey'
    )

    ez_dft = sim0.add_dft_fields(
        center=adj_source_center,
        size=adj_source_size,
        frequency=[fcen],
        field_component='Ez'
    )

    # register stop condition dft fields
    sim0.add_dft_fields(
        center=geom_center,
        size=geom_size,
        frequency=[fcen],
        field_component='Ex'
    )

    sim0.run(
        stop_condition=stop_condition,
        np=20,
        skip=True,
        pop=True
    )

    #%%
    adj_z, adj_y, adj_x = ffip.getgrid(adj_source_center, adj_source_size, adj_source_dim)

    ex_vals = ex_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ey_vals = ey_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ez_vals = ez_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    e_vals = np.sqrt(np.abs(ex_vals)**2 + np.abs(ey_vals)**2 + np.abs(ez_vals)**2)

    adj_dx = adj_x[1] - adj_x[0]
    adj_dy = adj_y[1] - adj_y[0]
    extent = (adj_x[0] - adj_dx/2, adj_x[-1] - adj_dx/2, adj_y[0] - adj_dy/2, adj_y[-1] - adj_dy/2)
    plt.imshow(np.squeeze(e_vals), origin='lower', extent=extent)
    plt.colorbar()
    plt.show()

    return e_vals, adj_x, adj_y, adj_z

e_vals_target, adj_x, adj_y, adj_z = get_target_fields()
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

adj_src = ffip.Adjoint_Source(
    adjoint_simulation=sim_adjoint,
    forward_simulation=sim_forward,
    function=src_func,
    frequency=fcen,
    center=adj_source_center,
    size=adj_source_size,
    dim=adj_source_dim,
    functionals=[
        ('|E|', e_vals_target)
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

scale = 1
rho_list = []
f_list = []
e_list = []

def fun(rho):
    adj_vol.density = np.ones(adj_vol.shape) * np.reshape(rho / scale, adj_vol.shape[1:])
    print('running forward calculation')
    sim_forward.run(
        stop_condition=stop_condition, 
        np=20,
        stdout=subprocess.DEVNULL
        # skip=True,
        # pop=True
    )
    
    return adj_src.eval_functionals_and_set_sources()

def fprime(rho):
    print('running adjoint calculation')
    sim_adjoint.run(
        stop_condition=stop_condition,
        np=20,
        stdout=subprocess.DEVNULL
        # skip=True,
        # pop=True
    )

    se = np.sum(adj_vol.get_sensitivity() / scale, axis=0)
    return se.ravel()

def show(rho):
    rho_t = np.reshape(rho / scale, adj_vol.shape[1:])

    rho_list.append(np.array(rho_t))

    f_list.append(adj_src.eval_functionals_and_set_sources())

    e_list.append(
        np.sqrt(
            np.abs(adj_src.forward_fields['Ex'])**2 + 
            np.abs(adj_src.forward_fields['Ey'])**2 + 
            np.abs(adj_src.forward_fields['Ez'])**2
            )
    )

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
        options={'maxiter' : 10, 'disp' : True}
    )

    print(res)

    with h5py.File(prefix + 'result.h5', 'w') as file:
        # save e rho
        file.create_dataset('rho target', data=bowtie_density)
        file.attrs['rho target x'] = bowtie_x
        file.attrs['rho target y'] = bowtie_y

        # save target e
        file.create_dataset('e target', data=e_vals_target)
        file.attrs['e target x'] = adj_src.x
        file.attrs['e target y'] = adj_src.y

        # save funs
        file.create_dataset('fun', data=np.array(f_list))

        # save rhos and es
        for i in range(len(rho_list)):
             file.create_dataset('rho %d' % i, data=rho_list[i])
             file.create_dataset('e %d' % i, data=e_list[i])
        

sensitivity_test()
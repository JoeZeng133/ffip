import ffip
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import minimize, Bounds

#%% setting up geometry

def disk(r=30):
    def _disk(x):
        return 1.0 * ((x[..., 1]**2 + x[..., 2]**2) < r**2)
    
    return _disk

def box(a=30, b=30):
    def _box(x):
        return 1.0 * (np.bitwise_and(np.abs(x[..., 1]) < b / 2, np.abs(x[..., 2]) < a / 2))
    return _box

#%% run target geometry
prefix = 'disk_'
dx = 4
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(200, 200, 100) + dpml * 2
fsrc = 1 / 500
fcen = 1 / 500

m1 = ffip.Au
m2 = ffip.Medium()

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

adj_source_size = ffip.Vector3(x=150, y=150)
adj_source_center = ffip.Vector3(z=30)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))

geom_size = ffip.Vector3(150, 150, 60)
geom_dim = (geom_size/dx+1).round()
geom_center = ffip.Vector3()

geom = ffip.Two_Medium_Box(
    size=geom_size,
    center=geom_center,
    dim=geom_dim,
    density=disk(r=50),
    medium1=m1,
    medium2=m2
)

geometry = [geom]

stop_condition = ffip.run_until_dft(
    center=geom_center,
    size=geom_size,
    field_component='Ex',
    time_interval_examined=4/fsrc,
    var=5e-3,
    frequency=fcen
)

sources = [ffip.Source(
    function=src_func,
    center=ffip.Vector3(z=-sim_size.z/2+dpml),
    size=ffip.Vector3(sim_size.x, sim_size.y, 0),
    field_component='Ex')
]

def get_target_fields():
    sim0 = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        pmls=pmls,
        geometry=geometry,
        sources=sources,
        progress_interval=20,
        input_file=prefix + 'target_input.h5',
        fields_output_file=prefix + 'target_output.h5'
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

    #%%
    sim0.run(
        stop_condition=stop_condition,
        np=2,
        skip=True,
        pop=True
    )

    #%%
    adj_z, adj_y, adj_x = ffip.getgrid(adj_source_center, adj_source_size, adj_source_dim)

    ex_vals = ex_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ey_vals = ey_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ez_vals = ez_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    e_vals = np.sqrt(np.abs(ex_vals)**2 + np.abs(ey_vals)**2 + np.abs(ez_vals)**2)
    
    plt.subplot(121)
    plt.imshow(np.squeeze(e_vals), origin='lower')
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(geom.density[0, ...], origin='lower')
    plt.show()

    return e_vals

e_vals = get_target_fields()

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
    input_file=prefix + 'adjoint_input.h5',
    fields_output_file=prefix + 'adjoint_output.h5',
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
        ('|E|', e_vals)
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

# rho in [0, scale] to control step size
scale = 1

def fun(rho):
    # broadcasting to build extrusion
    adj_vol.density = np.ones(adj_vol.shape) * np.reshape(rho / scale, adj_vol.shape[1:])
    print('running forward calculation')
    sim_forward.run(
        stop_condition=stop_condition, 
        np=2
        # skip=True,
        # pop=True
    )
    
    return adj_src.eval_functionals_and_set_sources()

def fprime(rho):
    print('running adjoint calculation')
    sim_adjoint.run(
        stop_condition=stop_condition,
        np=2
        # skip=True,
        # pop=True
    )

    se = np.sum(adj_vol.get_sensitivity() / scale, axis=0)
    return se.ravel()

def show(rho):
    rho_t = np.reshape(rho / scale, adj_vol.shape[1:])

    with h5py.File('bowtie_rho.h5', 'a') as file:
        tmp = len(file.keys())
        file.create_dataset('rho %d' % tmp, data=rho_t)
        file.close()
    
    plt.figure(1)
    plt.imshow(rho_t, origin='lower')
    plt.pause(3)

#%% sensitivity test
def sensitivity_test():
    rho0 = np.zeros(adj_vol.shape[1:]).flatten() + 0.5
    f1 = fun(rho0)
    print('f1=', f1)
    print('exp(f2 - f1)=', np.sum(fprime(rho0) * 0.01))
    sim_forward.fields_output_file = 'bowtie_forward_output1.h5'
    sim_forward.input_file = 'bowtie_forward_input1.h5'
    f2 = fun(rho0 + 0.01)
    print('f2=', f2)
    print('f2-f1=', f2 - f1)
    print("end")

#%% optimization
def optimize():
    rho0 = np.zeros(adj_vol.shape[1:]).flatten() + 0.5
    res = minimize(
        fun,
        rho0,
        method='L-BFGS-B',
        jac=fprime,
        bounds=Bounds(0, scale),
        callback=show,
        options={'maxiter' : 10, 'disp' : True}
    )

def level_set_optimize(fun, x0, jac, bounds):
    x_list = [x0]
    f_list = []
    fp_list = []
    itr = 0
    maxiter = 10
    # maximum step size
    ds = 0.5

    while 1:
        x = x_list[-1]
        # get functions and gradients
        f = fun(x)
        fp = np.reshape(jac(x), adj_vol.shape[1:])

        f_list.append(f)
        fp_list.append(fp)

        # get gradient of rho
        g0 = np.gradient(x, axis=0)
        g1 = np.gradient(x, axis=1)
        g = np.sqrt(g0**2 + g1**2)

        dx = g * fp
        tmp = np.max(np.abs(fp.ravel()))
        if tmp > ds:
            dx = ds / tmp  * dx
        
        x -= dx
        # clamp back to bounds
        x[x > bounds[1]] = bounds[1]
        x[x < bounds[0]] = bounds[0]

        x_list.append(x)
        plt.imshow(x, vmin=bounds[0], vmax=bounds[1], origin='lower')
        plt.show()
        print('at iter=%d, fun=%f' % (itr, f))

        itr += 1

        if itr > maxiter:
            break
    
    with h5py.File(prefix+'result.h5', 'w') as file:
        file.create_dataset('fun', data=np.array(f_list))
        for i in range(itr):
            file.create_dataset('x %d', data=x_list[i])
            file.create_dataset('xp %d', data=x_list[i])

pts = np.stack(np.meshgrid(adj_vol.z, adj_vol.y, adj_vol.x, indexing='ij'), axis=-1)
x0 = box(a=100, b=100)(pts)[0, ...]
plt.imshow(x0)
plt.show()

level_set_optimize(
    fun=fun,
    x0=x0,
    jac=fprime,
    bounds=(0, 1)
)

        






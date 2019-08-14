import ffip
import matplotlib.pyplot as plt
import numpy as np
import nlopt
import h5py
import subprocess
from filter_gen import get_green_filter, get_qT, q_inc, get_norm_green_filter
from initial_guess import get_rho0, get_bowtie_den
from scipy.signal import convolve

#%% basic setup
prefix = 'rectangular_temp_'
dx = 4
e_inc = 103
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(500, 500, 100) + dpml * 2
fsrc = 1 / 500
fcen = 1 / 800

#materials and susceptibilti to build permittivity
m1 = ffip.Au
m2 = ffip.Medium()
fept = ffip.FePt(fcen)

sus1 = ffip.Au_susc[0]
sus2 = ffip.Au_susc[3]

e1 = m1.get_epsilon(fcen)
e2 = m2.get_epsilon(fcen)

#nonlinear permittivity interpolation
def epsilon_fun(rho):
    return (np.sqrt(e1) * rho + np.sqrt(e2) * (1 - rho))**2

def epsilon_der(rho):
    return 2 * (np.sqrt(e1) * rho + np.sqrt(e2) * (1 - rho)) * (np.sqrt(e1) - np.sqrt(e2))

# get source input DTFT for use in normalization
src_func = ffip.Gaussian1(fsrc, start_time=0.5/fcen)
ref_t = np.arange(src_func.start_time, src_func.end_time, dt)
ref_fft = np.sum(np.exp(-2j*np.pi*fcen*ref_t) * src_func(ref_t))

pmls = [ffip.PML(
    thickness=dpml
)]

# plane wave source comming from below
sources = [ffip.Source(
    function=src_func,
    center=ffip.Vector3(z=-sim_size.z/2+dpml),
    size=ffip.Vector3(sim_size.x, sim_size.y, 0),
    field_component='Ex')
]

# adjoint source region (objective function influence region)
adj_source_size = ffip.Vector3(x=300, y=300)
adj_source_center = ffip.Vector3(z=25)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))
adj_z, adj_y, adj_x = ffip.getgrid(center=adj_source_center, size=adj_source_size, dim=adj_source_dim)

green_filter = get_norm_green_filter(
    q_inc, e_inc,
    dx=adj_x[1] - adj_x[0],
    dy=adj_y[1] - adj_y[0],
    lt=1, 
    dimt=100
)

# adjoint volume region where density can vary
geom_size = ffip.Vector3(200, 200, 60)
geom_dim = (geom_size/dx+1).round()
geom_shape = (int(geom_dim.z), int(geom_dim.y), int(geom_dim.x))
geom_center = ffip.Vector3(z=-10)
geom_z, geom_y, geom_x = ffip.getgrid(center=geom_center, size=geom_size, dim=geom_dim)


bowtie = ffip.General_Medium_Box(
    size=geom_size,
    center=geom_center,
    dim=geom_dim,
    density=get_bowtie_den(),
    epsilon_fun=epsilon_fun,
    frequency=fcen,
    e_sus=[sus1, sus2],
    suffix='0'
)

# stop until dft of the adjoint region converges
stop_condition = ffip.run_until_dft(
    center=geom_center,
    size=geom_size,
    field_component='Ex',
    time_interval_examined=4/fsrc,
    var=5e-3,
    frequency=fcen
)

# unperturbed geometry
fept_layer = ffip.Box(
    size=ffip.Vector3(x=sim_size.x, y=sim_size.y, z=20),
    center=ffip.Vector3(z=30),
    material=fept
)

gold_layer = ffip.Box(
    size=ffip.Vector3(x=sim_size.x, y=sim_size.y, z=geom_size.z),
    center=geom_center,
    material=m1
)

def get_incident_fields():
    sim0 = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        pmls=pmls,
        sources=sources,
        # geometry=[gold_layer, fept_layer],
        input_file=prefix + "incident_input.h5",
        fields_output_file=prefix+"incident_output.h5",
        progress_interval=20
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
        skip=True,
        # pop=True,
        stop_condition=stop_condition,
        np=35
    )

    adj_z, adj_y, adj_x = ffip.getgrid(adj_source_center, adj_source_size, adj_source_dim)

    ex_vals = ex_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ey_vals = ey_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ez_vals = ez_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    e_vals = np.sqrt(np.abs(ex_vals)**2 + np.abs(ey_vals)**2 + np.abs(ez_vals)**2)

    T = convolve(np.squeeze(e_vals)**2, green_filter, mode='same')
    plt.subplot(121)
    plt.imshow(np.squeeze(e_vals), extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
    plt.subplot(122)
    plt.imshow(T, origin='lower', extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
    plt.show()

    return e_vals, T

def get_max_fields():
    sim0 = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        pmls=pmls,
        sources=sources,
        geometry=[bowtie, gold_layer, fept_layer],
        input_file=prefix + "max_input.h5",
        fields_output_file=prefix+"max_output.h5",
        progress_interval=20
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
        skip=True,
        # pop=True,
        stop_condition=stop_condition,
        np=35
    )

    adj_z, adj_y, adj_x = ffip.getgrid(adj_source_center, adj_source_size, adj_source_dim)

    ex_vals = ex_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ey_vals = ey_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ez_vals = ez_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    e_vals = np.sqrt(np.abs(ex_vals)**2 + np.abs(ey_vals)**2 + np.abs(ez_vals)**2)

    T = convolve(np.squeeze(e_vals)**2, green_filter, mode='same')
    print('emax=', np.max(e_vals)/e_inc, ',emin=', np.min(e_vals)/e_inc)
    
    plt.subplot(121)
    plt.imshow(np.squeeze(e_vals) / e_inc, extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(T, origin='lower', extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
    plt.colorbar()
    plt.show()

    return e_vals, T

e_vals, T = get_max_fields()
# T_incident = T[T.shape[0]//2, T.shape[1]//2]

#%% adjoint setups

# weighted sum of magnitdue
def udf1_weight(weight):

    def obj1(ex, ey, ez):
        e = np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2)

        res = np.sum(weight * e)

        return -res

    def fpconj1(ex, ey, ez):
        e = np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2)

        pex = -np.squeeze(ex) / e * weight
        pey = -np.squeeze(ey) / e * weight
        pez = -np.squeeze(ez) / e * weight

        return (np.reshape(np.conj(pex), ex.shape),
            np.reshape(np.conj(pey), ey.shape),
            np.reshape(np.conj(pez), ez.shape))
    
    return (['Ex', 'Ey', 'Ez'], fpconj1, obj1)

def udf2_e_match(e_target):

    def obj2(ex, ey, ez):
        e = np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2)

        res = 0.5 * np.sum((e - e_target)**2)

        return res

    def fpconj2(ex, ey, ez):
        e = np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2)

        dfde = e - e_target

        pex = ex / e * dfde
        pey = ey / e * dfde
        pez = ez / e * dfde

        return (np.reshape(np.conj(pex), ex.shape),
            np.reshape(np.conj(pey), ey.shape),
            np.reshape(np.conj(pez), ez.shape))
    
    return (['Ex', 'Ey', 'Ez'], fpconj2, obj2)

def udf3_t_match(t_target, mask):

    def obj3(ex, ey, ez):
        q = np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2
        t = convolve(q, green_filter, 'same')
        res = 0.5 * np.sum(((t - t_target)*mask)**2)

        return res

    def fpconj3(ex, ey, ez):
        e = np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2)
        t = convolve(e**2, green_filter, 'same')

        dfdt = (t - t_target)*mask
        dfdq = convolve(dfdt, np.flip(green_filter), 'same')
        dfde = 2 * e * dfdq

        pex = ex / e * dfde
        pey = ey / e * dfde
        pez = ez / e * dfde

        return (np.reshape(np.conj(pex), ex.shape),
            np.reshape(np.conj(pey), ey.shape),
            np.reshape(np.conj(pez), ez.shape))
    
    return (['Ex', 'Ey', 'Ez'], fpconj3, obj3)

def get_gaussian_filter(r = 3, sigma=1.5):
    filter_x = np.arange(-r, r+1)
    filter_X, filter_Y = np.meshgrid(filter_x, filter_x, indexing='ij')
    dfilter = np.exp(-0.5 * (filter_X**2 + filter_Y**2) / sigma)
    return dfilter

def filt(x, dfilter):
    return ffip.TO_convolve(x, dfilter)

def filt_transpose(xp, dfilter):
    return ffip.TO_convolve_transpose(xp, dfilter)

def thresh(x, beta, eta):
    return (np.tanh(beta * eta) + np.tanh(beta * (x - eta))) / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))

def thresh_transpose(xp, x, beta, eta):
    return xp * beta * (1 / np.cosh(beta * (x - eta)))**2 / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))

rho_list = []
rho_filt_list = []
rho_thresh_list = []
f_list = []
e_list = []
t_list = []
se_list = []

def get_fun(udf, dfilter, beta=1, alpha=0.01, nsc=5, eta=0.5, save=False):

    # objective function and derivatives, target temperature mutliple times of incident temperature
    sim_forward = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        sources=sources,
        geometry=[fept_layer, gold_layer],
        pmls=pmls,
        input_file=prefix + 'forward_input.h5',
        fields_output_file=prefix + 'forward_output.h5',
        progress_interval=20
    )

    sim_adjoint = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        geometry=[fept_layer, gold_layer],
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
        udfs=[udf]
    )

    adj_vol = ffip.Adjoint_Volume(
        adjoint_simulation=sim_adjoint,
        forward_simulation=sim_forward,
        frequency=fcen,
        center=geom_center,
        size=geom_size,
        dim=geom_dim,
        density=None,
        norm=ref_fft,
        epsilon_fun=epsilon_fun,
        epsilon_der=epsilon_der,
        e_sus=[sus1, sus2]
    )

    iter = 0
    cur_f = np.inf
    acc_nsc = 0

    def fun(x, grad):
        nonlocal iter
        nonlocal cur_f
        nonlocal acc_nsc

        # standard method
        rho = np.reshape(x, geom_shape[1:])
        rho_filt = filt(rho, dfilter)
        rho_thresh = thresh(rho_filt, beta=beta, eta=eta)

        adj_vol.density = np.ones(adj_vol.shape) * rho_thresh
        print('running forward calculation')
        sim_forward.run(
            stdout=subprocess.DEVNULL,
            stop_condition=stop_condition, 
            np=35
        )

        res = adj_src.eval_functionals_and_set_sources()
        
        print('running adjoint calculation')
        sim_adjoint.run(
            stdout=subprocess.DEVNULL,
            stop_condition=stop_condition,
            np=35
        )

        if grad.size > 0:
            se = np.sum(adj_vol.get_sensitivity(), axis=0)
            se_thresh = thresh_transpose(se, rho_filt, beta=beta, eta=eta)
            se_filt = filt_transpose(se_thresh, dfilter)

            grad[:] = se_filt.ravel()

        # save result by identifying iterations as improving evaluation
        if save and res < cur_f:

            iter += 1
            f_list.append(res)
            e_list.append(
                np.sqrt(
                    np.abs(adj_src.forward_fields['Ex'])**2 + 
                    np.abs(adj_src.forward_fields['Ey'])**2 + 
                    np.abs(adj_src.forward_fields['Ez'])**2
                    )
            )

            t = convolve(np.squeeze(e_list[-1]**2), green_filter, mode='same')

            t_list.append(t)

            rho_list.append(np.array(rho))
            rho_filt_list.append(np.array(rho_filt))
            rho_thresh_list.append(np.array(rho_thresh))

            se_list.append(se_filt)

            plt.figure(1)
            plt.clf()

            plt.subplot(131)
            plt.imshow(rho_thresh, origin='lower', vmin=0, vmax=1, extent=(geom_x[0], geom_x[-1], geom_y[0], geom_y[-1]))
            plt.colorbar()
            plt.subplot(132)
            plt.imshow(t, origin='lower', extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
            plt.colorbar()
            plt.subplot(133)
            plt.imshow(np.squeeze(e_list[-1]), origin='lower', extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
            plt.colorbar()

            plt.pause(3)
        
            print("at iter", iter, ",fun=", res)

            if np.abs((cur_f - res) / res) < alpha:
                acc_nsc += 1
                print("Solution doesnt change much after %d iterations" % acc_nsc)
                if acc_nsc >= nsc:
                    raise ValueError("Solution doesnt change much after %d iterations" % acc_nsc)
            
            else:
                acc_nsc = 0
            
            cur_f = res

        return res
    
    return fun

#%% optimization 1 run
def optimize(rho0, filename, udf, dfilter, beta=1, alpha=0.02, nsc=100):

    global rho_list
    global rho_filt_list
    global rho_thresh_list
    global f_list
    global e_list
    global t_list
    global se_list

    rho_list = []
    rho_filt_list = []
    rho_thresh_list = []
    f_list = []
    e_list = []
    t_list = []
    se_list = []

    try:
        opt = nlopt.opt(nlopt.LD_MMA, rho0.size)
        opt.set_lower_bounds(0)
        opt.set_upper_bounds(1)
        opt.set_min_objective(get_fun(udf, dfilter=dfilter, beta=beta, alpha=alpha, nsc=nsc, eta=0.5, save=True))
        opt.set_maxeval(100)
        opt.optimize(rho0.flatten())    

    except ValueError:
        pass

    with h5py.File(filename, 'w') as file:
        # save metadata
        file.create_dataset('rho x', data=geom_x)
        file.create_dataset('rho y', data=geom_y)
        file.create_dataset('t x', data=adj_x)
        file.create_dataset('t y', data=adj_y)

        for i in range(len(rho_list)):
            file.create_dataset('rho %d' % i, data=rho_list[i])
            file.create_dataset('rho thresh %d' % i, data=rho_thresh_list[i])
            file.create_dataset('rho filt %d' % i, data=rho_filt_list[i])
            file.create_dataset('e %d' % i, data=e_list[i])
            file.create_dataset('se %d' % i, data=se_list[i])
            file.create_dataset('t %d' % i, data=t_list[i])

        file.create_dataset('fun', data=np.array(f_list))

#%% setting up
rho0 = np.ones(geom_shape[1:])
r = 2
beta = 1
dfilter = get_gaussian_filter(3, r / 2)
filename = prefix + 'result.h5'
qnorm, T = get_qT(lx=adj_source_size.x, ly=adj_source_size.y, lt=1, dimx=adj_source_dim.x, dimy=adj_source_dim.y, dimt=100)
e_target = np.sqrt(qnorm) * e_inc

plt.imshow(e_target, origin='lower', extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
plt.show()
# udf = udf1_weight(np.reshape(q, adj_source_shape))
udf = udf2_e_match(np.reshape(e_target, adj_source_shape))

#%% sensitivity test
def sensitivity_test():
    fun = get_fun(udf=udf, dfilter=get_gaussian_filter(), save=False)

    rho0 = np.zeros(geom_shape[1:]) + 0.5
    print('rho shape=', rho0.shape)

    grad = np.zeros(rho0.size)
    f1 = fun(rho0, grad)

    print('f1=', f1)
    print('exp(f2 - f1)=', np.sum(grad * 0.005))

    f2 = fun(rho0 + 0.005, grad)
    print('f2=', f2)
    print('f2-f1=', f2 - f1)
    print("end")


#%%
# sensitivity_test()

# with h5py.File(prefix + 'result0.h5', 'r') as file:
#     numiter = float(file.attrs['numiter'])
#     rho0 = np.array(file['rho %d' % (numiter - 1)])
#     beta = file.attrs['beta'] * 1.8

for i in range(0, 5):
    filename = prefix+'result%d.h5' % i
    print('beta=', beta, ',r=', r, '#####################')
    optimize(rho0=rho0, filename=filename, udf=udf, dfilter=dfilter, beta=beta, alpha=0.01, nsc=3)

    with h5py.File(filename, 'a') as file:
        file.attrs['beta'] = beta
        file.attrs['numiter'] = len(rho_list)
    
    beta = beta * 1.8
    rho0 = np.array(rho_list[-1])

plt.show()
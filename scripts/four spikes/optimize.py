import ffip
import matplotlib.pyplot as plt
import numpy as np
import nlopt
import h5py
import subprocess
from filter_gen import get_green_filter
from scipy.signal import convolve

#%% basic setup
prefix = 'rectangular_temp_'
dx = 4
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(500, 500, 160) + dpml * 2
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

adj_source_size = ffip.Vector3(x=200, y=200)
adj_source_center = ffip.Vector3(z=30)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))
green_filter = get_green_filter(
    lx=adj_source_size.x, 
    ly=adj_source_size.y, 
    lt=1, 
    dimx=adj_source_dim.x, 
    dimy=adj_source_dim.y, 
    dimt=100
)

geom_size = ffip.Vector3(200, 200, 60)
geom_dim = (geom_size/dx+1).round()
geom_center = ffip.Vector3(z=-10)

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
    size=ffip.Vector3(x=sim_size.x, y=sim_size.y, z=60),
    center=ffip.Vector3(z=-10),
    material=m1
)

def get_incident_fields():
    sim0 = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        pmls=pmls,
        sources=sources,
        geometry=[fept_layer, gold_layer],
        input_file=prefix + "incident_input.h5",
        fields_output_file=prefix+"incident_output.h5",
        progress_interval=4
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
        # skip=True,
        # pop=True,
        stop_condition=stop_condition,
        np=30
    )

    adj_z, adj_y, adj_x = ffip.getgrid(adj_source_center, adj_source_size, adj_source_dim)

    ex_vals = ex_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ey_vals = ey_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ez_vals = ez_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    e_vals = np.sqrt(np.abs(ex_vals)**2 + np.abs(ey_vals)**2 + np.abs(ez_vals)**2)

    T = convolve(np.squeeze(e_vals)**2, green_filter, mode='same')
    plt.imshow(T, origin='lower', extent=(adj_x[0], adj_x[-1], adj_y[0], adj_y[-1]))
    plt.show()

    return e_vals, T, adj_x, adj_y

e_vals, T, adj_x, adj_y = get_incident_fields()
#%% adjoint setups
spotx = 30
spoty = 60
py, px = np.meshgrid(adj_y, adj_x, indexing='ij')
mask = (np.abs(px) < spotx / 2) and (np.abs(py) < spoty / 2)
T_center = T[T.shape[0]//2, T.shape[1]//2]

def objective(ex, ey, ez, mask = None, T_target=0):
    shape = ex.shape

    if mask is None:
        mask = np.ones(shape[1:])

    e2 = np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2

    t = convolve(np.squeeze(e2), green_filter, mode='same')

    res = 0.5 * np.sum((t - T_target)**2 * mask)

    return res


def fex_conj(ex, ey, ez, mask = None, T_target=0):
    if mask is None:
        mask = np.ones(ex.shape[1:])

    e = np.squeeze(np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2))

    t = convolve(e**2, green_filter, mode='same')

    #df/dt
    pt = (t - T_target) * mask
    #df/dq
    pq = convolve(pt, np.flip(green_filter), mode='same')
    #df/de
    pe = 2 * e * pq
    #df/dex
    pew = np.squeeze(ex) / e * pe

    return np.reshape(np.conj(pew), ex.shape)

def fey_conj(ex, ey, ez, mask = None, T_target=0):
    if mask is None:
        mask = np.ones(ex.shape[1:])

    e = np.squeeze(np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2))

    t = convolve(e**2, green_filter, mode='same')

    #df/dt
    pt = (t - T_target) * mask
    #df/dq
    pq = convolve(pt, np.flip(green_filter), mode='same')
    #df/de
    pe = 2 * e * pq
    #df/dex
    pew = np.squeeze(ex) / e * pe

    return np.reshape(np.conj(pew), ex.shape)

def fez_conj(ex, ey, ez, mask = None, T_target=0):
    if mask is None:
        mask = np.ones(ex.shape[1:])

    e = np.squeeze(np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2))

    t = convolve(e**2, green_filter, mode='same')

    #df/dt
    pt = (t - T_target) * mask
    #df/dq
    pq = convolve(pt, np.flip(green_filter), mode='same')
    #df/de
    pe = 2 * e * pq
    #df/dex
    pew = np.squeeze(ex) / e * pe

    return np.reshape(np.conj(pew), ex.shape)


# sim_forward = ffip.Simulation(
#     size=sim_size,
#     resolution=1/dx,
#     sources=sources,
#     geometry=[fept_layer, gold_layer],
#     pmls=pmls,
#     input_file=prefix + 'forward_input.h5',
#     fields_output_file=prefix + 'forward_output.h5',
#     progress_interval=20
# )

# sim_adjoint = ffip.Simulation(
#     size=sim_size,
#     resolution=1/dx,
#     pmls=pmls,
#     input_file=prefix+'adjoint_input.h5',
#     fields_output_file=prefix+'adjoint_output.h5',
#     progress_interval=20
# )

# adj_src = ffip.Adjoint_Source(
#     adjoint_simulation=sim_adjoint,
#     forward_simulation=sim_forward,
#     function=src_func,
#     frequency=fcen,
#     center=adj_source_center,
#     size=adj_source_size,
#     dim=adj_source_dim,
#     udfs=(['Ex', 'Ey', 'Ez'], [fex_conj, fey_conj, fez_conj), objective])
# )

# adj_vol = ffip.Adjoint_Volume(
#     adjoint_simulation=sim_adjoint,
#     forward_simulation=sim_forward,
#     frequency=fcen,
#     center=geom_center,
#     size=geom_size,
#     dim=geom_dim,
#     density=None,
#     # medium1=m1,
#     # medium2=m2,
#     norm=ref_fft,
#     epsilon_fun=epsilon_fun,
#     epsilon_der=epsilon_der,
#     e_sus=[sus1, sus2]
# )

# def get_gaussian_filter(r = 5):

#     filter_sigma = r / 2
#     filter_x = np.arange(-r, r+1)
#     filter_X, filter_Y = np.meshgrid(filter_x, filter_x, indexing='ij')
#     filter = np.exp(-0.5 * (filter_X**2 + filter_Y**2) / filter_sigma)
#     filter = filter / np.sum(filter.ravel())

#     return filter

# def filt(x, filter):
#     return convolve(x, filter, mode='valid')

# def filt_transpose(xp, filter):
#     tp = np.flip(filter, axis=None)
#     return convolve(xp, tp, mode='full')

# def thresh(x, beta, eta):
#     return (np.tanh(beta * eta) + np.tanh(beta * (x - eta))) / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))

# def thresh_transpose(xp, x, beta, eta):
#     return xp * beta * (1 / np.cosh(beta * (x - eta)))**2 / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))

# rho_list = []
# rho_filt_list = []
# rho_thresh_list = []
# f_list = []
# e_list = []
# se_list = []

# filter = get_gaussian_filter(3)
# rho0 = np.zeros(adj_vol.shape[1:]) + 0.5
# rho0 = np.pad(rho0, (3, 3), 'edge')

# # with h5py.File(prefix + 'previous.h5', 'r') as file:
# #     rho0 = np.array(file['rho 44'])

# with h5py.File(prefix + 'result.h5', 'w') as file:
#     # save e rho
#     file.create_dataset('rho target', data=rho_target)
#     file.create_dataset('rho x', data=bowtie_x)
#     file.create_dataset('rho y', data=bowtie_y)

#     # save target e
#     file.create_dataset('e target', data=e_vals_target)
#     file.create_dataset('e x', data=adj_x)
#     file.create_dataset('e y', data=adj_y)


# def get_fun(filter, beta=1, alpha=0.01, nsc=10, eta=0.5, maxiter = 15, save=False):
#     iter = 0
#     cur_f = np.inf

#     def fun(x, grad):
#         nonlocal iter
#         nonlocal cur_f

#         # standard method
#         rho = np.reshape(x, rho0.shape)
#         rho_filt = filt(rho, filter)
#         rho_thresh = thresh(rho_filt, beta=beta, eta=eta)

#         adj_vol.density = np.ones(adj_vol.shape) * rho_thresh
#         print('running forward calculation')
#         sim_forward.run(
#             stdout=subprocess.DEVNULL,
#             stop_condition=stop_condition, 
#             np=3
#         )

#         res = adj_src.eval_functionals_and_set_sources()
        
#         print('running adjoint calculation')
#         sim_adjoint.run(
#             stdout=subprocess.DEVNULL,
#             stop_condition=stop_condition,
#             np=3
#         )

#         if grad.size > 0:
#             se = np.sum(adj_vol.get_sensitivity(), axis=0)
#             se_thresh = thresh_transpose(se, rho_filt, beta=beta, eta=eta)
#             se_filt = filt_transpose(se_thresh, filter)

#             grad[:] = se_filt.ravel()

#         # save result by identifying iterations as improving evaluation
#         if save and res < cur_f:

#             iter += 1
#             f_list.append(res)
#             e_list.append(
#                 np.sqrt(
#                     np.abs(adj_src.forward_fields['Ex'])**2 + 
#                     np.abs(adj_src.forward_fields['Ey'])**2 + 
#                     np.abs(adj_src.forward_fields['Ez'])**2
#                     )
#             )

#             rho_list.append(np.array(rho))
#             rho_filt_list.append(np.array(rho_filt))
#             rho_thresh_list.append(np.array(rho_thresh))

#             se_list.append(se_filt)

#             plt.figure(1)
#             plt.imshow(rho_thresh, origin='lower', vmin=0, vmax=1)
#             plt.pause(3)
        
#             print("at iter", iter, ",fun=", res)

#             if iter > nsc and np.abs(cur_f - res) / np.abs(res) < alpha:
#                 raise ValueError("Beta Increase Needed")
            
#             cur_f = res

#         return res
    
#     return fun
    
# #%% sensitivity test
# def sensitivity_test():
#     fun = get_fun(filter=get_gaussian_filter(r=3), beta=3)
#     rho0 = np.zeros(adj_vol.shape[1:]) + 0.5
#     rho0 = np.pad(rho0, (3, 3), 'edge')

#     grad = np.zeros(rho0.size)
#     f1 = fun(rho0, grad)

#     print('f1=', f1)
#     print('exp(f2 - f1)=', np.sum(grad * 0.003))
#     sim_forward.fields_output_file = prefix+'forward_pert_output.h5'
#     sim_forward.input_file = prefix+'forward_pert_input.h5'

#     f2 = fun(rho0 + 0.003, grad)
#     print('f2=', f2)
#     print('f2-f1=', f2 - f1)
#     print("end")

# #%%
# def optimize():
#     beta = 1

#     for i in range(8):
#         lastiter = len(rho_list)
#         try:
#             print("#################### current beta=", beta)
#             opt = nlopt.opt(nlopt.LD_MMA, rho0.size)
#             opt.set_lower_bounds(0)
#             opt.set_upper_bounds(1)
#             opt.set_min_objective(get_fun(filter=filter, beta=beta, alpha=0.02, nsc=10, eta=0.5, save=True))
#             opt.set_maxeval(30)
#             if len(f_list) > 0:
#                 opt.optimize(rho_list[-1].flatten())
#             else:
#                 opt.optimize(rho0.flatten())    

#         except ValueError:
#             pass
        
#         beta = beta * 1.5

#         with h5py.File(prefix + 'result.h5', 'a') as file:
#             # save rhos and es
#             for i in range(lastiter, len(rho_list)):
#                 file.create_dataset('rho %d' % i, data=rho_list[i])
#                 file.create_dataset('rho thresh %d' % i, data=rho_thresh_list[i])
#                 file.create_dataset('rho filt %d' % i, data=rho_filt_list[i])
#                 file.create_dataset('e %d' % i, data=e_list[i])
#                 file.create_dataset('se %d' % i, data=se_list[i])
        
#         lastiter = len(f_list)

#     with h5py.File(prefix + 'result.h5', 'a') as file:
#         # save funs
#         file.create_dataset('fun', data=np.array(f_list))

# sensitivity_test()
# # optimize()
# plt.show()
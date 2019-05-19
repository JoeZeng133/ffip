import ffip
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, Bounds

#%% setting up geometry

def get_bowtie_den():
    a = 30
    b = 40
    c = 50
    d = 50

    tri1 = ffip.Polygon([(0, b), (c, b + d), (-c, b + d)])
    tri2 = ffip.Polygon([(0, -b), (c, -b-d), (-c, -b-d)])
    rect = ffip.box(-a/2, -b-d, a/2, b+d)

    res = tri1.union(tri2).union(rect)

    x, dx = np.linspace(-100, 100, 50, retstep=True)
    y, dy = np.linspace(-100, 100, 50, retstep=True)
    z = 0
    pts = np.stack(np.meshgrid(z, y, x, indexing='ij'), axis=-1)

    density = ffip.planar_polygon(res)
    rho = density(pts)

    extent = (x[0] - dx/2, x[-1] - dx/2, y[0] - dy/2, y[-1] - dy/2)
    plt.imshow(np.squeeze(rho), origin='lower', extent=extent)
    plt.colorbar()
    plt.xlabel('x [nm]')
    plt.ylabel('y [nm]')
    plt.show()
    return density

density = get_bowtie_den()

#%% run target geometry
dx = 2
dt = 0.5 * dx
dpml = 8 * dx
sim_size = ffip.Vector3(300, 300, 100) + dpml * 2
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

adj_source_size = ffip.Vector3(x=100, y=100)
adj_source_center = ffip.Vector3(z=40)
adj_source_dim = (adj_source_size/dx+1).round()
adj_source_shape = (int(adj_source_dim.z), int(adj_source_dim.y), int(adj_source_dim.x))

geom_size = ffip.Vector3(200, 200, 60)
geom_dim = (geom_size/dx+1).round()
geom_center = ffip.Vector3()

bowtie = ffip.Two_Medium_Box(
    size=geom_size,
    center=geom_center,
    dim=geom_dim,
    density=density,
    medium1=m1,
    medium2=m2
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

sources = [ffip.Source(
    function=src_func,
    center=ffip.Vector3(z=-sim_size.z/2+dpml),
    size=ffip.Vector3(sim_size.x, sim_size.y, 0),
    field_component='Ex')
]

def get_bowtie_fields():
    sim0 = ffip.Simulation(
        size=sim_size,
        resolution=1/dx,
        pmls=pmls,
        geometry=geometry,
        sources=sources,
        progress_interval=20,
        input_file='bowtie_input.h5',
        fields_output_file='bowtie_output.h5'
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
        skip=True)

    #%%
    adj_z, adj_y, adj_x = ffip.getgrid(adj_source_center, adj_source_size, adj_source_dim)

    ex_vals = ex_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ey_vals = ey_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)
    ez_vals = ez_dft(fcen, adj_z, adj_y, adj_x).reshape(adj_source_shape)

    return ex_vals, ey_vals, ez_vals

ex_vals, ey_vals, ez_vals = get_bowtie_fields()
e_vals = np.sqrt(np.abs(ex_vals)**2 + np.abs(ey_vals)**2 + np.abs(ez_vals)**2)

plt.imshow(np.squeeze(e_vals), origin='lower')
plt.colorbar()
plt.show()

#%% adjoint setups
sim_forward = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    sources=sources,
    pmls=pmls,
    input_file='bowtie_forward_input.h5',
    fields_output_file='bowtie_forward_output.h5',
    progress_interval=20
)

sim_adjoint = ffip.Simulation(
    size=sim_size,
    resolution=1/dx,
    pmls=pmls,
    input_file='bowtie_adjoint_input.h5',
    fields_output_file='bowtie_adjoint_output.h5',
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

rho0 = np.zeros(adj_vol.shape[1:])
scale = 1

def func(rho):
    adj_vol.density = np.ones(adj_vol.shape) * np.reshape(rho, adj_vol.shape[1:]) / scale
    print('running forward calculation')
    sim_forward.run(
        stop_condition=stop_condition, 
        np=2,
        skip=True,
        pop=True
    )
    
    return adj_src.eval_functionals_and_set_sources()

def fprime(rho):
    print('running adjoint calculation')
    sim_adjoint.run(
        stop_condition=stop_condition,
        np=2,
        skip=True,
        pop=True
    )

    se = np.sum(adj_vol.get_sensitivity(), axis=0) / scale
    return se.ravel()

def show(rho):
    se = np.sum(adj_vol.get_sensitivity(), axis=0)
    plt.imshow(se, origin='lower')
    plt.show()

#%%
minimize(
    func,
    rho0,
    method='L-BFGS-B',
    jac=fprime,
    bounds=Bounds(0, scale),
    callback=show
)
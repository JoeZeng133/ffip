import ffip
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import minimize, Bounds

#%% setting up geometry

def get_bowtie_den():
    a = 30
    b = 60
    c = 400
    d = 400

    tri1 = ffip.Polygon([(a/2, b/2), (c/2, d/2), (c/2, -d/2), (a/2, -b/2)])
    tri2 = ffip.Polygon([(-a/2, b/2), (-c/2, d/2), (-c/2, -d/2), (-a/2, -b/2)])

    res = tri1.union(tri2)

    x, dx = np.linspace(-300, 300, 121, retstep=True)
    y, dy = np.linspace(-300, 300, 121, retstep=True)
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
    return density, np.squeeze(rho), extent

def read_heat_gen():
    file = h5py.File('heat_gen.h5')
    x = np.array(file['x'])[:, 0]
    y = np.array(file['y'])[0, :]
    q = np.array(file['q'])
    T = np.array(file['T'])

    dx = x[1] - x[0]
    dy = y[1] - y[0]
    extent = (x[0] - dx/2, x[-1] - dx/2, y[0] - dy/2, y[-1] - dy/2)

    plt.subplot(121)
    plt.imshow(q, origin='lower', extent=extent)
    plt.xlabel('x [nm]')
    plt.ylabel('y [nm]')
    plt.subplot(122)
    plt.imshow(T, origin='lower', extent=extent)
    plt.xlabel('x [nm]')
    plt.ylabel('y [nm]')
    plt.show()

    return x, y, q, T

density_fun, rho_bowtie, rho_bowtie_extent = get_bowtie_den()
qx, qy, q, T = read_heat_gen()

#%% run target geometry
prefix = 'four_spikes_'
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

adj_source_size = ffip.Vector3(x=100, y=100)
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

adj_dx = adj_x[1] - adj_x[0]
adj_dy = adj_y[1] - adj_y[0]
extent = (adj_x[0] - adj_dx/2, adj_x[-1] - adj_dx/2, adj_y[0] - adj_dy/2, adj_y[-1] - adj_dy/2)
plt.imshow(np.squeeze(e_vals), origin='lower', extent=extent)
plt.colorbar()
plt.show()
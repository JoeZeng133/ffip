import ffip
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import minimize, Bounds
import subprocess

def get_bowtie_den():
    a = 10
    b = 10
    c = 220
    d = 220

    tri1 = ffip.Polygon([(a/2, b/2), (c/2, d/2), (c/2, -d/2), (a/2, -b/2)])
    tri2 = ffip.Polygon([(-a/2, b/2), (-c/2, d/2), (-c/2, -d/2), (-a/2, -b/2)])

    res = tri1.union(tri2)
    density = ffip.planar_polygon(res)
    # pts = np.stack(np.meshgrid(0, y, x, indexing='ij'), axis=-1)

    return density

def get_rec_aperture_den(y, x):
    a=30
    b=60
    res = ffip.Polygon([(-a/2, b/2), (a/2, b/2), (a/2, -b/2), (-a/2, -b/2)])

    density = ffip.planar_polygon(res)
    pts = np.stack(np.meshgrid(0, y, x, indexing='ij'), axis=-1)

    return 1 - density(pts)

def get_rho0(y, x):
    return np.squeeze(get_rec_aperture_den(y, x)) * 1.0

if __name__ == "__main__":
    dx = 4
    geom_size = ffip.Vector3(200, 200, 60)
    geom_dim = (geom_size/dx+1).round()
    geom_shape = (int(geom_dim.z), int(geom_dim.y), int(geom_dim.x))
    geom_center = ffip.Vector3(z=-10)
    geom_z, geom_y, geom_x = ffip.getgrid(center=geom_center, size=geom_size, dim=geom_dim)

    rho0 = get_rho0(geom_y, geom_x)

    print(rho0.shape)

    plt.imshow(rho0, origin='lower', extent=(geom_x[0], geom_x[-1], geom_y[0], geom_y[-1]))
    plt.show()
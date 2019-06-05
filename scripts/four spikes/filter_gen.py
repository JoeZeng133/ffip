import sympy
import numpy as np
import matplotlib.pyplot as plt    
import h5py
from scipy.signal import convolve
from scipy.optimize import minimize, Bounds

# time scale, lenght scale, rescaled diffusivity
t = 1e-9
l = 1e-9
a = 4.6440e-07 / l**2 * t

def green(x, t):
    return 1 / (4 * np.pi * a * t) * np.exp(-np.sum(x**2, axis=-1) / (4 * a * t))

def get_green_filter(lx=200, ly=200, lt=1, dimx=51, dimy=51, dimt=100):

    dx = lx / (dimx - 1)
    dy = ly / (dimy - 1)
    dt = lt / (dimt - 1)

    x = np.linspace(-lx/2, lx/2, dimx)
    y = np.linspace(-ly/2, ly/2, dimy)
    pts = np.stack(np.meshgrid(y, x, indexing='ij'), axis=-1)

    filter = np.zeros((int(dimy), int(dimx)))
    filter[int(dimy//2), int(dimx//2)] = dt

    for i in range(1, dimt):
        filter += green(pts, i * dt) * dx * dy * dt

    return filter

if __name__ == "__main__":
    lx = 200
    ly = 200
    lt = 1
    dimx = 51
    dimy = 51
    dimt = 100
    dx = lx / (dimx - 1)
    dy = ly / (dimy - 1)
    dt = lt / (dimt - 1)

    x = np.linspace(-lx/2, lx/2, dimx)
    y = np.linspace(-ly/2, ly/2, dimy)
    pts = np.stack(np.meshgrid(y, x, indexing='ij'), axis=-1)

    extent = (-lx/2 - dx/2, lx/2 - dx/2, -ly/2 - dy/2, ly/2 - dy/2)
    filter = get_green_filter(lx=lx, ly=ly, lt=lt, dimx=dimx, dimy=dimy, dimt=dimt)

    def rectangular(x, ax=30, ay=60):
        return 1.0 * (np.abs(x[..., 0]) < ay/2) * (np.abs(x[..., 1]) < ax/2)

    mask = rectangular(pts, ax=10, ay=20)
    num = np.sum(mask)
    T_target = rectangular(pts, ax=30, ay=60)

    def fun(x):
        tmp = np.reshape(x, T_target.shape)
        T = convolve(tmp, filter, mode='same')
        Tmask = T * mask
        ave = np.sum(Tmask) / num

        res =  0.5 * np.sum(((Tmask - ave) * mask)**2, axis=None)
        return res

    def fp(x):
        tmp = np.reshape(x, T_target.shape)
        T = convolve(tmp, filter, mode='same')
        Tmask = T * mask
        ave = np.sum(Tmask) / num

        # df/dT
        res = ((Tmask - ave) - 1 / num * np.sum((T - ave)*mask)) * mask
        # df/dq
        res = convolve(res, np.flip(filter), mode='same')

        return res.flatten()

    res = minimize(
        fun=fun,
        x0=np.ones(T_target.size),
        jac=fp,
        method='L-BFGS-B',
        bounds=Bounds(1, np.inf)
    )

    print(res)

    q = np.reshape(res.x, T_target.shape)
    T = convolve(q, filter, mode='same')
    plt.imshow(T, origin='lower')
    plt.show()

    with h5py.File('filter_result.h5', 'w') as file:
        file.create_dataset('q', data=q)
        file.create_dataset('T', data=T)
        file.create_dataset('T target', data=T_target)
        file.create_dataset('filter', data=filter)
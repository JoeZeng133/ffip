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
k = 1.5 * l * t
q_scale = l**3 * t
q_inc = 8e17 * q_scale * 5

def green_fun(x, t):
    return 1 / (4 * np.pi * a * t) * np.exp(-np.sum(x**2, axis=-1) / (4 * a * t))

def get_green_filter(dx=4, dy=4, lt=1, dimt=100):

    # truncate the gaussian function at 3 sigma
    sigma = np.sqrt(2 * a * lt)
    R = 3 * sigma

    dimx = round(R / dx) * 2 + 1
    lx = (dimx - 1) * dx / 2
    dt = lt / (dimt - 1)

    x = np.linspace(-lx/2, lx/2, dimx)
    pts = np.stack(np.meshgrid(x, x, indexing='ij'), axis=-1)

    filter = np.zeros((int(dimx), int(dimx)))
    filter[int(dimx//2), int(dimx//2)] = dt

    for i in range(1, dimt):
        filter += green_fun(pts, i * dt) * dx * dy * dt

    return filter

def get_norm_green_filter(q_inc, e_inc, dx=4, dy=4, lt=1, dimt=100):
    gfilter = get_green_filter(dx=dx, dy=dy, lt=lt, dimt=dimt)

    return a / k * gfilter * q_inc / e_inc**2

def rectangular_profile(x, ax=30, ay=60):
    return 1.0 * (np.abs(x[..., 0]) < ay/2) * (np.abs(x[..., 1]) < ax/2)

def get_qT(lx=200, ly=200, lt=1, dimx=51, dimy=51, dimt=100):

    dx = lx / (dimx - 1)
    dy = ly / (dimy - 1)

    x = np.linspace(-lx/2, lx/2, dimx)
    y = np.linspace(-ly/2, ly/2, dimy)
    pts = np.stack(np.meshgrid(y, x, indexing='ij'), axis=-1)

    filter = get_green_filter(dx=dx, dy=dy, lt=lt, dimt=dimt)

    w = 1e-3

    T_mask = rectangular_profile(pts, ax=30, ay=60)
    T_target =  T_mask * 500
    # mask = rectangular_profile(pts, ax=60, ay=60)
    mask = T_mask

    def fun(x):
        tmp = np.reshape(x, T_target.shape)
        T = convolve(tmp, filter, mode='same')

        res =  0.5 * np.sum(((T - T_target) * mask)**2, axis=None) + w * np.sum((tmp * (1 - mask))**2)
        return res

    def fp(x):
        tmp = np.reshape(x, T_target.shape)
        T = convolve(tmp, filter, mode='same')

        # df/dT
        res = (T - T_target) * mask
        # df/dq
        res = convolve(res, np.flip(filter), mode='same')

        return res.flatten() + w * 2 * x * (1 - mask.flatten())
    
    s_inc = q_inc * (a/k)
    lb = s_inc * ((8.3/103)**2)
    ub = s_inc * (2.55**2)

    print('s_inc=', s_inc, ',s_lb=', lb, ',s_ub=', ub)
 
    res = minimize(
        fun=fun,
        x0=np.ones(T_target.size) * lb,
        jac=fp,
        method='L-BFGS-B',
        bounds=Bounds(lb, ub)
    )

    print(res)

    s = np.reshape(res.x, T_target.shape)
    T = convolve(s, filter, mode='same')

    qphys = s / (a / k) / q_scale
    qnorm = s / s_inc

    print("maximum q=", np.max(qphys), ',maximum qnorm=', np.max(qnorm))

    plt.subplot(121)
    plt.imshow(T, origin='lower', extent=(x[0], x[-1], y[0], y[-1]))
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(np.sqrt(qnorm), origin='lower', extent=(x[0], x[-1], y[0], y[-1]))
    plt.colorbar()
    plt.show()

    return qnorm, T

if __name__ == "__main__":
    lx = 300
    ly = 300
    lt = 1
    dimx = 151
    dimy = 151
    dimt = 100

    get_qT(lx, ly, lt, dimx, dimy, dimt)
import ffip
import matplotlib.pyplot as plt
import numpy as np
import nlopt
import h5py
import subprocess
from scipy.signal import convolve
import numpy.random
import scipy.misc

def get_gaussian_filter(r = 3):

    filter_sigma = r / 2
    filter_x = np.arange(-r, r+1)
    filter_X, filter_Y = np.meshgrid(filter_x, filter_x, indexing='ij')
    filter = np.exp(-0.5 * (filter_X**2 + filter_Y**2) / filter_sigma)

    return filter

def to_filter_test():
    filter = np.random.rand(5, 5)
    # filter = get_gaussian_filter(r=5)
    p1 = np.ones((100, 100))
    dp = np.zeros((100, 100))
    dp[-1, -1] = 100

    pf1 = ffip.TO_convolve(p1, filter)
    f1 = np.sum(pf1)

    fp = np.ones((100, 100))
    fp1 = ffip.TO_convolve_transpose(fp, filter)
    print("exp diff=", np.sum(fp1 * dp))

    p2 = p1 + dp
    pf2 = ffip.TO_convolve(p2, filter)
    f2 = np.sum(pf2)

    print("act diff=", f2 - f1)

def np_filter_test():
    filter = np.random.rand(5, 5)
    p1 = np.ones((100, 100))
    dp = np.zeros((100, 100))
    dp[-1, -1] = 100

    pf1 = convolve(p1, filter, 'same')
    f1 = np.sum(pf1)

    fp = np.ones((100, 100))
    fp1 = convolve(fp, np.flip(filter), 'same')
    print("exp diff=", np.sum(fp1 * dp))

    p2 = p1 + dp
    pf2 = convolve(p2, filter, 'same')
    f2 = np.sum(pf2)

    print("act diff=", f2 - f1)

np_filter_test()
import ffip
import matplotlib.pyplot as plt
import numpy as np
import nlopt
import h5py
import subprocess
from scipy.signal import convolve


m1 = ffip.Au_susc
m2 = ffip.Au_LD_susc

for sus in m1:
    print("gamma/omega=", sus.gamma / sus.frequency, ",e=", sus.get_epsilon(1/800) / sus.sigma, ",sigma=", sus.sigma)

for sus in m2:
    print("gamma/omega=", sus.gamma/sus.frequency, ",e=", sus.get_epsilon(1/800) / sus.sigma, ",sigma=", sus.sigma)
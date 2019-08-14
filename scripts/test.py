import ffip
import matplotlib.pyplot as plt
import numpy as np
import nlopt
import h5py
import subprocess
from scipy.signal import convolve
import numpy.random
import scipy.misc


a = np.array([1,2,3,4])
print(a[a>2].size)
import matplotlib.pyplot as plt
import numpy as np


for i in range(3):
    plt.figure(i)
    plt.imshow(np.random.random((30, 30)))
    plt.draw()
    plt.pause(1)

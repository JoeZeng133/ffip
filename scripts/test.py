import numpy as np
import ffip
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
import matplotlib.pyplot as plt

x = np.linspace(-6, 6, 100)
y = np.linspace(-5, 5, 100)
z = 100
pts = np.stack(np.meshgrid(z, y, x, indexing='ij'), axis=-1)

p1 = Polygon(((-1.0, -1.0), (-1.0, 1.0), (1.0, 1.0), (1.0, -1.0)))
p2 = Polygon(((-1, 2), (1, 2), (0, 4)))

res = p1.union(p2)

den = ffip.planar_polygon(res, center=ffip.Vector3(x=0.5, z=100), height=1)
rho = den(pts)

extend = (x[0], x[-1], y[-1], y[1])
plt.imshow(np.squeeze(rho), origin='lower', extent=extend)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


import math
import numbers
from math import cos, sin, sqrt, ceil, floor
from numbers import Number
from copy import deepcopy
import numpy as np
from abc import ABCMeta, abstractmethod
import ffip
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, box

def cmp_shape(shape1=(), shape2=()):
    if len(shape1) != len(shape2):
        return False
    
    for i in range(len(shape1)):
        if shape1[i] != shape2[i]:
            return False
    
    return True

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

def metadata_to_pts(*args):
    # meshgrid and stack xn, xn-1, ..., x1 together
    return np.stack(np.meshgrid(*args, indexing='ij'), axis=-1)

class Vector3(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = float(x) if type(x) is int else x
        self.y = float(y) if type(y) is int else y
        self.z = float(z) if type(z) is int else z

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __add__(self, other):
        if (isinstance(other, Number)):
            other = Vector3(other, other, other)

        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z

        return Vector3(x, y, z)

    def __sub__(self, other):
        if (isinstance(other, Number)):
            other = Vector3(other, other, other)

        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z

        return Vector3(x, y, z)

    def __mul__(self, other):
        if type(other) is Vector3:
            return self.dot(other)
        elif isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError(
                "No operation known for 'Vector3 * {}'".format(type(other)))

    def __truediv__(self, other):
        if type(other) is Vector3:
            return Vector3(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, Number):
            return Vector3(self.x / other, self.y / other, self.z / other)
        else:
            raise TypeError(
                "No operation known for 'Vector3 / {}'".format(type(other)))

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError(
                "No operation known for '{} * Vector3'".format(type(other)))

    def __getitem__(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        else:
            raise IndexError("No value at index {}".format(i))

    def __setitem__(self, i, data):
        if i == 0:
            self.x = data
        elif i == 1:
            self.y = data
        elif i == 2:
            self.z = data
        else:
            raise IndexError("No value at index {}".format(i))

    def __repr__(self):
        return "Vector3<{}, {}, {}>".format(self.x, self.y, self.z)

    def __array__(self):
        return np.array([self.x, self.y, self.z])

    def conj(self):
        return Vector3(self.x.conjugate(), self.y.conjugate(), self.z.conjugate())

    def scale(self, s):
        x = self.x * s
        y = self.y * s
        z = self.z * s

        return Vector3(x, y, z)

    def dot(self, v):
        return self.x * v.x + self.y * v.y + self.z * v.z

    def cdot(self, v):
        return self.conj().dot(v)

    def cross(self, v):
        x = self.y * v.z - self.z * v.y
        y = self.z * v.x - self.x * v.z
        z = self.x * v.y - self.y * v.x

        return Vector3(x, y, z)

    def norm(self):
        return math.sqrt(abs(self.dot(self)))

    def unit(self):
        return self.scale(1 / self.norm())
    
    def prod(self):
        return self.x * self.y * self.z

    def close(self, v, tol=1.0e-7):
        return (abs(self.x - v.x) <= tol and
                abs(self.y - v.y) <= tol and
                abs(self.z - v.z) <= tol)

    def rotate(self, axis, theta):
        u = axis.unit()
        vpar = u.scale(u.dot(self))
        vcross = u.cross(self)
        vperp = self - vpar
        return vpar + (vperp.scale(math.cos(theta)) + vcross.scale(math.sin(theta)))

    def get_json(self):
        return [self.x, self.y, self.z]

    def ceil(self):
        return Vector3(ceil(self.x), ceil(self.y), ceil(self.z))

    def floor(self):
        return Vector3(floor(self.x), floor(self.y), floor(self.z))

    def round(self):
        return Vector3(round(self.x), round(self.y), round(self.z))
    
    def copy(self):
        return Vector3(self.x, self.y, self.z)


class Medium:

    def __init__(self,
                 epsilon=1,
                 mu=1,
                 E_conductivity=0,
                 M_conductivity=0,
                 E_susceptibilities=[],
                 M_susceptibilities=[]):

        self.epsilon = float(epsilon)
        self.mu = float(mu)
        self.E_susceptibilities = deepcopy(E_susceptibilities)
        self.M_susceptibilities = deepcopy(M_susceptibilities)
        self.E_conductivity = float(E_conductivity)
        self.M_conductivity = float(M_conductivity)

    # return wave impedance of non dispersive part
    def get_imp(self):
        return sqrt(self.mu / self.epsilon)

    # return wave speed of non dispersive part
    def get_c(self):
        return 1 / sqrt(self.mu * self.epsilon)

    def get_json(self):

        e_sus = [item.get_json() for item in self.E_susceptibilities]
        m_sus = [item.get_json() for item in self.M_susceptibilities]

        return {
            'epsilon': self.epsilon,
            'mu': self.mu,
            'electric conductivity': self.E_conductivity,
            'magnetic conductivity': self.M_conductivity,
            'electric susceptibility': e_sus,
            'magnetic susceptibility': m_sus
        }

    # return epsilon at a particular frequency
    def get_epsilon(self, frequency):
        res = self.epsilon - 1j * self.E_conductivity / (2 * np.pi * frequency)
        for item in self.E_susceptibilities:
            res += item.get_epsilon(frequency)

        return res
    
    def get_dis_epsilon(self, frequency, dt):
        z = np.exp(1j * 2 * np.pi * frequency * dt)
        res = self.epsilon + 2 * dt * z / (z**2 - 1) * self.E_conductivity

        for item in self.E_susceptibilities:
            res += item.get_dis_epsilon(frequency, dt)
        
        return res

class Susceptibility:

    def __init__(self, sigma=1.0):
        self.sigma = sigma

    def build_eqs(self, a0=0.0, a1=0.0, b0=0.0, b1=0.0, b2=0.0):
        self.a0 = float(a0)
        self.a1 = float(a1)
        self.b0 = float(b0)
        self.b1 = float(b1)
        self.b2 = float(b2)

    def __eq__(self, other):
        return (self.a0 == other.a0 and
                self.a1 == other.a1 and
                self.b0 == other.b0 and
                self.b1 == other.b1 and
                self.b2 == other.b2)

    def get_epsilon(self, frequency):
        w = 2 * np.pi * frequency
        return self.sigma * (self.a0 + self.a1 * 1j * w) / (self.b0 + 1j * w * self.b1 - w * w * self.b2)

    def get_dis_epsilon(self, frequency, dt):
        c0 = 2 * self.b2 + self.b1 * dt
        c1 = (4 * self.b2 - 2 * self.b0 * dt * dt) / c0
        c2 = (-2 * self.b2 + self.b1 * dt) / c0
        c3 = (2 * self.a0 * dt * dt) / c0

        z = np.exp(1j * 2 * np.pi * frequency * dt)

        return self.sigma * c3 * z / (z**2 - c1 * z - c2)




class LorentzianSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        w = 2 * np.pi * frequency
        g = 2 * np.pi * gamma
        super().__init__(**kwargs)
        super().build_eqs(a0=w**2, b0=w**2, b1=g, b2=1)
        self.frequency = float(frequency)
        self.gamma = float(gamma)

    def get_json(self):
        return {"type": "Lorentz",
                "frequency": self.frequency,
                "gamma": self.gamma,
                "amplitude": self.sigma}


class DrudeSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        w = 2 * np.pi * frequency
        g = 2 * np.pi * gamma
        super().__init__(**kwargs)
        super().build_eqs(a0=w**2, b1=g, b2=1)
        self.frequency = float(frequency)
        self.gamma = float(gamma)

    def get_json(self):
        return {"type": "Drude",
                "frequency": self.frequency,
                "gamma": self.gamma,
                "amplitude": self.sigma}

class DeybeSusceptibility(Susceptibility):

    def __init__(self, tau=0.0, **kwargs):
        super().__init__(**kwargs)
        super().build_eqs(a0=1, b0=1, b1=tau)
        self.tau = float(tau)

    
    def get_json(self):
        return {'type' : 'Deybe',
                'tau' : self.tau,
                'amplitude' : self.sigma
                }

class GeometricObject:

    def __init__(self, material=Medium(), center=Vector3()):
        self.material = deepcopy(material)
        self.center = center.copy()

    def __add__(self, vec):
        self.center += vec

    def shift(self, vec):
        c = deepcopy(self)
        c.center += vec
        return c


class Sphere(GeometricObject):

    def __init__(self, radius, **kwargs):
        self.radius = float(radius)
        super().__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, val):
        self._radius = float(val) if val > 0 else 0.0
    
    def dump(self, file):
        pass

    def get_json(self):
        return {'type': 'sphere',
                'center': self.center.get_json(),
                'radius': self.radius,
                'medium': self.material.get_json()}


class Box(GeometricObject):

    def __init__(self, size=Vector3(), **kwargs):
        self.size = size
        super().__init__(**kwargs)

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, val):
        self._size = val.copy() if (val.x > 0 and val.y >
                             0 and val.z > 0) else Vector3()
    
    def dump(self, file):
        pass

    def get_json(self):
        return {"type": "box",
                "center": self.center.get_json(),
                "size": self.size.get_json(),
                "medium": self.material.get_json()}

class Two_Medium_Box:
    def __init__(self, size=Vector3(), center=Vector3(), dim=Vector3(), density=None, medium1=Medium(), medium2=Medium(), suffix='0'):
        # suffix adjoint is reserved
        self.size = size.copy()
        self.center = center.copy()
        self.dim = dim.round()
        self.medium1 = deepcopy(medium1)
        self.medium2 = deepcopy(medium2)

        if callable(density):
            self.density = density(np.stack(np.meshgrid(self.z, self.y, self.x, indexing='ij'), axis=-1))

        elif isinstance(density, np.ndarray):
            self.density = density

        elif density is None:
            self.density = np.zeros(self.shape, dtype=float)
        else:
            raise ValueError('amplitude input is not numpy array or a function')
        
        self.dataset_name = "two medium box %s" % suffix
    
    @property
    def density(self):
        return self._density
    
    @density.setter
    def density(self, val):
        if not cmp_shape(self.shape, val.shape):
            raise ValueError("density size is not compatible with the box dimension")

        self._density = np.array(val, copy=True)
    
    @property
    def x(self):
        return np.linspace(self.center.x - self.size.x/2, self.center.x + self.size.x/2, self.dim.x)
    
    @property
    def y(self):
        return np.linspace(self.center.y - self.size.y/2, self.center.y + self.size.y/2, self.dim.y)

    @property
    def z(self):
        return np.linspace(self.center.z - self.size.z/2, self.center.z + self.size.z/2, self.dim.z)

    @property
    def dimension(self):
        return self.dim
    
    @property
    def shape(self):
        return (int(self.dim.z), int(self.dim.y), int(self.dim.x))

    @property
    def numel(self):
        return int(self.dim.prod())
    
    def dump(self, file):
        file.create_dataset(self.dataset_name, data=self.density.ravel())

    def get_json(self):
        return {"type" : "two medium box",
                "center" : self.center.get_json(),
                "size" : self.size.get_json(),
                "dimension" : self.dim.get_json(),
                "density dataset" : self.dataset_name,
                "medium1" : self.medium1.get_json(),
                "medium2" : self.medium2.get_json()
                }

def getgrid(center=Vector3(), size=Vector3(), dim=Vector3()):
    
    x = np.linspace(center.x - size.x/2, center.x + size.x/2, dim.x)
    y = np.linspace(center.y - size.y/2, center.y + size.y/2, dim.y)
    z = np.linspace(center.z - size.z/2, center.z + size.z/2, dim.z)

    return z, y, x

def planar_polygon(polygon, center=Vector3(), height=np.inf):

    def one_geom(g, req_pts):
        xy_check = ffip.check_inside(
                req_pts=np.stack((req_pts[..., 2] - center.x, req_pts[..., 1] - center.y), axis=-1), 
                poly_pts=np.asarray(g.exterior)
            )
        z_check = np.abs(req_pts[..., 0] - center.z) < (height/2)

        return xy_check * z_check

    def internal(req_pts):
        res = False
        if isinstance(polygon, MultiPolygon):
            for g in polygon.geoms:
                res = np.bitwise_or(res, one_geom(g, req_pts)) 
        else:
            res = one_geom(polygon, req_pts)
    
        return res * 1

    return internal

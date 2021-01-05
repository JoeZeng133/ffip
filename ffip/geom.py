import math
import numbers
from math import cos, sin, sqrt, ceil, floor
from numbers import Number
from copy import deepcopy
import numpy as np
from numpy.linalg import inv
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


class Susceptibility:

    def __init__(self, sigma=1.0):
        self.sigma = sigma

    def build_eqs(self, b0=0.0, b1=0.0, a0=0.0, a1=0.0, a2=0.0):
        self.b0 = float(b0)
        self.b1 = float(b1)
        self.a0 = float(a0)
        self.a1 = float(a1)
        self.a2 = float(a2)

    def __eq__(self, other):
        return (self.b0 == other.b0 and
                self.b1 == other.b1 and
                self.a0 == other.a0 and
                self.a1 == other.a1 and
                self.a2 == other.a2)

    def get_epsilon(self, frequency):
        s = 2j * np.pi * frequency
        return self.sigma * (self.b0 + self.b1 * s) / (self.a0 + self.a1 * s + self.a2 * s**2)

    def get_dis_epsilon(self, frequency, dt):
        a0 = self.a0 + 2 * self.a1 / dt + 4 * self.a2 / dt / dt
        a1 = 2 * self.a0 - 8 * self.a2 / dt / dt
        a2 = self.a0 - 2 * self.a1 / dt + 4 * self.a2 / dt / dt
        b0 = self.b0 + 2 * self.b1 / dt
        b1 = 2 * self.b0
        b2 = self.b0 - 2 * self.b1 / dt

        z = np.exp(1j * 2 * np.pi * frequency * dt)

        return self.sigma * (b0 + b1 * z**-1 + b2 * z**-2) / (a0 + a1 * z**-1 + a2 * z**-2)


class IIR_Susceptibility(Susceptibility):

    def __init__(self, a: complex, c: complex, **kwargs):
        self.a = complex(a)
        self.c = complex(c)
        super().__init__(**kwargs)
        super().build_eqs(
            -2 * np.real(a * np.conj(c)),
            2 * np.real(c),
            abs(a)**2,
            -2 * np.real(a),
            1)

    def get_json(self):
        return {"type": "IIR",
                "b0": self.b0,
                "b1": self.b1,
                "a0": self.a0,
                "a1": self.a1,
                "a2": self.a2,
                "amplitude": self.sigma}


class LorentzianSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        w = 2 * np.pi * frequency
        g = 2 * np.pi * gamma
        super().__init__(**kwargs)
        super().build_eqs(b0=w**2, a0=w**2, a1=g, a2=1)
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
        super().build_eqs(b0=w**2, a1=g, a2=1)
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
        super().build_eqs(b0=1, a0=1, a1=tau)
        self.tau = float(tau)

    def get_json(self):
        return {'type': 'Deybe',
                'tau': self.tau,
                'amplitude': self.sigma
                }


class Medium:

    def __init__(self,
                 epsilon: float = 1,
                 mu: float = 1,
                 E_conductivity: float = 0,
                 M_conductivity: float = 0,
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
            self.density = density(np.stack(np.meshgrid(
                self.z, self.y, self.x, indexing='ij'), axis=-1))

        elif isinstance(density, np.ndarray):
            self.density = density

        elif density is None:
            self.density = np.zeros(self.shape, dtype=float)
        else:
            raise ValueError(
                'amplitude input is not numpy array or a function')

        self.dataset_name = "two medium box %s" % suffix

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, val):
        if not cmp_shape(self.shape, val.shape):
            raise ValueError(
                "density size is not compatible with the box dimension")

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
        return {"type": "two medium box",
                "center": self.center.get_json(),
                "size": self.size.get_json(),
                "dimension": self.dim.get_json(),
                "density dataset": self.dataset_name,
                "medium1": self.medium1.get_json(),
                "medium2": self.medium2.get_json()
                }


class Param_Medium:
    # parameterized medium
    def __init__(
        self,
        medium: Medium,
        e_inf_fun,
        esus_amp_fun
    ):
        self.medium = deepcopy(medium)
        self.e_inf_fun = deepcopy(e_inf_fun)
        self.esus_amp_fun = deepcopy(esus_amp_fun)
        # force only e_inf and e_sus
        self.medium.mu = float(1)
        self.medium.E_conductivity = float(0)
        self.medium.M_conductivity = float(0)
        self.medium.M_susceptibilities = []
        # enable easy calcuation of susceptibility with unit amplitude
        for sus in self.medium.E_susceptibilities:
            sus.sigma = 1.0

    def get_epsilon(self, rho: np.ndarray, frq):
        sus = [sus.get_epsilon(frq)
                for sus in self.medium.E_susceptibilities]
        return self.e_inf(rho) + np.sum(self.esus_amp(rho) * np.array(sus), axis=-1)

    def get_dis_epsilon(self, rho: np.ndarray, frq, dt):
        sus = [sus.get_dis_epsilon(frq, dt)
                for sus in self.medium.E_susceptibilities]
        return self.e_inf(rho) + np.sum(self.esus_amp(rho) * np.array(sus), axis=-1)

    def get_medium(self, rho: float) -> Medium:
        esus_amp = self.esus_amp_fun(rho)
        self.medium.epsilon = float(self.e_inf_fun(rho))

        for i in range(len(self.medium.E_susceptibilities)):
            self.medium.E_susceptibilities[i].sigma = esus_amp[i]

        return deepcopy(self.medium)

    def e_inf(self, rho: np.ndarray) -> np.ndarray:
        return self.e_inf_fun(rho)

    def m_inf(self, rho: np.ndarray) -> np.ndarray:
        return np.ones(rho.shape)

    def esus_amp(self, rho: np.ndarray) -> np.ndarray:
        return self.esus_amp_fun(rho)

    def msus_amp(self, rho: np.ndarray) -> np.ndarray:
        return np.zeros(0)


class Param_Medium_Box:

    def __init__(
            self,
            size=Vector3,
            center=Vector3,
            density=None,
            dim=Vector3(),
            medium: Param_Medium = None,
            suffix="0"):

        self.size = size.copy()
        self.center = center.copy()
        self.dim = dim.round()
        self.medium = deepcopy(medium)
        # sigma is included in get_epsilon routine which needs to be taken out

        if callable(density):
            self.density = density(np.stack(np.meshgrid(
                self.z, self.y, self.x, indexing='ij'), axis=-1))

        elif isinstance(density, np.ndarray):
            self.density = density

        elif density is None:
            self.density = np.zeros(self.shape, dtype=float)
        else:
            raise ValueError(
                'amplitude input is not numpy array or a function')

        prefix = "general medium box2 %s " % suffix

        self.density_dataset = prefix + "density"
        self.density_fun_dataset = prefix + "density function"

        self.epsilon_dataset = prefix + "epsilon"
        self.mu_dataset = prefix + "mu"
        self.e_sus_amp_dataset = prefix + "e sus amp"
        self.m_sus_amp_dataset = prefix + "m sus amp"

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, val):
        if not cmp_shape(self.shape, val.shape):
            raise ValueError(
                "density size is not compatible with the box dimension")

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
        # create interpolation sequence
        rho = np.linspace(0, 1, 5001)

        file.create_dataset(self.density_dataset, data=self.density.ravel())
        file.create_dataset(self.density_fun_dataset, data=rho.ravel())

        # einf and minf
        file.create_dataset(self.epsilon_dataset,
                            data=self.medium.e_inf(rho).ravel())
        file.create_dataset(
            self.mu_dataset, data=self.medium.m_inf(rho).ravel())

        # esus and msus
        file.create_dataset(self.e_sus_amp_dataset,
                            data=self.medium.esus_amp(rho).ravel())
        file.create_dataset(self.m_sus_amp_dataset,
                            data=self.medium.msus_amp(rho).ravel())

    def get_json(self):
        return {"type": "general medium box",
                "center": self.center.get_json(),
                "size": self.size.get_json(),
                "dimension": self.dim.get_json(),
                "density dataset": self.density_dataset,
                "density function dataset": self.density_fun_dataset,
                "epsilon dataset": self.epsilon_dataset,
                "mu dataset": self.mu_dataset,
                "electric susceptibility amplitudes dataset": self.e_sus_amp_dataset,
                "magnetic susceptibility amplitudes dataset": self.m_sus_amp_dataset,
                "medium": self.medium.get_json()
                }


class General_Medium_Box:

    def __init__(self, size=Vector3(), center=Vector3(), dim=Vector3(), frequency=1.0, density=None, epsilon_fun=None, e_sus=[], suffix='0'):

        self.size = size.copy()
        self.center = center.copy()
        self.dim = dim.round()
        self.frequency = float(frequency)

        if len(e_sus) != 2:
            raise ValueError("Currently only support two poles")

        self.e_sus = deepcopy(e_sus)
        self.medium = Medium(epsilon=1, E_susceptibilities=self.e_sus)
        self.epsilon_fun = deepcopy(epsilon_fun)
        # sigma is included in get_epsilon routine which needs to be taken out
        self.e1 = self.e_sus[0].get_epsilon(frequency) / self.e_sus[0].sigma
        self.e2 = self.e_sus[1].get_epsilon(frequency) / self.e_sus[1].sigma

        self.mat = np.array(
            [[np.real(self.e1), np.real(self.e2)],
             [np.imag(self.e1), np.imag(self.e2)]]
        )

        self.invmat = inv(self.mat)

        if callable(density):
            self.density = density(np.stack(np.meshgrid(
                self.z, self.y, self.x, indexing='ij'), axis=-1))

        elif isinstance(density, np.ndarray):
            self.density = density

        elif density is None:
            self.density = np.zeros(self.shape, dtype=float)
        else:
            raise ValueError(
                'amplitude input is not numpy array or a function')

        prefix = "general medium box %s " % suffix

        self.density_dataset = prefix + "density"
        self.density_fun_dataset = prefix + "density function"

        self.epsilon_dataset = prefix + "epsilon"
        self.mu_dataset = prefix + "mu"
        self.e_sus_amp_dataset = prefix + "e sus amp"
        self.m_sus_amp_dataset = prefix + "m sus amp"

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, val):
        if not cmp_shape(self.shape, val.shape):
            raise ValueError(
                "density size is not compatible with the box dimension")

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
        # create interpolation sequence
        rho = np.linspace(0, 1, 1000)
        epsilon = self.epsilon_fun(rho)

        epsilon_arr = np.stack(
            (np.real(epsilon) - 1, np.imag(epsilon)), axis=0)
        rho_arr = self.invmat @ epsilon_arr
        rho_arr = np.transpose(rho_arr)

        file.create_dataset(self.density_dataset, data=self.density.ravel())
        file.create_dataset(self.density_fun_dataset, data=rho.ravel())

        # einf = muinf = 1
        file.create_dataset(self.epsilon_dataset, data=np.ones(rho.size))
        file.create_dataset(self.mu_dataset, data=np.ones(rho.size))

        file.create_dataset(self.e_sus_amp_dataset, data=rho_arr.flatten())
        file.create_dataset(self.m_sus_amp_dataset, data=np.ones(0))

    def get_json(self):
        return {"type": "general medium box",
                "center": self.center.get_json(),
                "size": self.size.get_json(),
                "dimension": self.dim.get_json(),
                "density dataset": self.density_dataset,
                "density function dataset": self.density_fun_dataset,
                "epsilon dataset": self.epsilon_dataset,
                "mu dataset": self.mu_dataset,
                "electric susceptibility amplitudes dataset": self.e_sus_amp_dataset,
                "magnetic susceptibility amplitudes dataset": self.m_sus_amp_dataset,
                "medium": self.medium.get_json()
                }


def getgrid(center=Vector3(), size=Vector3(), dim=Vector3()):

    x = np.linspace(center.x - size.x/2, center.x + size.x/2, dim.x)
    y = np.linspace(center.y - size.y/2, center.y + size.y/2, dim.y)
    z = np.linspace(center.z - size.z/2, center.z + size.z/2, dim.z)

    return z, y, x


def planar_polygon(polygon, center=Vector3(), height=np.inf):

    def one_geom(g, req_pts):
        xy_check = ffip.check_inside(
            req_pts=np.stack(
                (req_pts[..., 2] - center.x, req_pts[..., 1] - center.y), axis=-1),
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

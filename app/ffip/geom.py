import math
import numbers
from math import cos, sin, sqrt, ceil, floor
from numbers import Number
from copy import deepcopy
import numpy as np
from abc import ABCMeta, abstractmethod

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

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
            raise TypeError("No operation known for 'Vector3 * {}'".format(type(other)))

    def __truediv__(self, other):
        if type(other) is Vector3:
            return Vector3(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, Number):
            return Vector3(self.x / other, self.y / other, self.z / other)
        else:
            raise TypeError("No operation known for 'Vector3 / {}'".format(type(other)))

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError("No operation known for '{} * Vector3'".format(type(other)))

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

class Medium:

    def __init__(self,
                 epsilon=1,
                 mu=1,
                 E_conductivity=0,
                 M_conductivity=0,
                 E_susceptibilities=[],
                 M_susceptibilities=[]):

        self.epsilon = epsilon
        self.mu = mu
        self.E_susceptibilities = E_susceptibilities
        self.M_susceptibilities = M_susceptibilities
        self.E_conductivity = E_conductivity
        self.M_conductivity = M_conductivity
    
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
            'epsilon' : self.epsilon, 
            'mu' : self.mu, 
            'electric conductivity' : self.E_conductivity, 
            'magnetic conductivity' : self.M_conductivity,
            'electric susceptibility' : e_sus,
            'magnetic susceptibility' : m_sus
            }
    
    # return epsilon at a particular frequency
    def get_epsilon(self, frequency):
        res = self.epsilon - 1j * self.E_conductivity / (2 * np.pi * frequency)
        for item in self.E_susceptibilities:
            res += item.get_epsilon(frequency)
        
        return res

class Susceptibility:

    def __init__(self, sigma=1.0):
        self.sigma = sigma

    def build_eqs(self, a0=0.0, a1=0.0, b0=0.0, b1=0.0, b2=0.0):
        self.a0 = a0
        self.a1 = a1
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2

    def __eq__(self, other):
        return (self.a0 == other.a0 and
                self.a1 == other.a1 and
                self.b0 == other.b0 and
                self.b1 == other.b1 and
                self.b2 == other.b2)
    
    def get_epsilon(self, frequency):
        w = 2 * np.pi * frequency
        return self.sigma * (self.a0 + self.a1 * 1j * w) / (self.b0 + 1j * w * self.b1 - w * w * self.b2)


class LorentzianSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        w = 2 * np.pi * frequency
        g = 2 * np.pi * gamma
        super().__init__(**kwargs)
        super().build_eqs(a0=w**2, b0=w**2, b1=g, b2=1)
        self.frequency = frequency
        self.gamma = gamma
    
    def get_json(self):
        return {"type" : "Lorentz", 
                "frequency" : self.frequency, 
                "gamma" : self.gamma, 
                "amplitude" : self.sigma}

class DrudeSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        w = 2 * np.pi * frequency
        g = 2 * np.pi * gamma
        super().__init__(**kwargs)
        super().build_eqs(a0=w**2, b1=g, b2=1)
        self.frequency = frequency
        self.gamma = gamma
    
    def get_json(self):
        return {"type" : "Drude", 
                "frequency" : self.frequency, 
                "gamma" : self.gamma, 
                "amplitude" : self.sigma}

class CriticalPointsSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, phi=0.0, **kwargs):
        super().__init__(**kwargs)
        super().build_eqs(
                        a0=2*frequency(frequency*cos(phi) - gamma*sin(phi)),
                        a1=-2*frequency*sin(phi),
                        b0=frequency*frequency+gamma*gamma,
                        b1=2*gamma,
                        b2=1)

class GeometricObject(object):

    def __init__(self, material=Medium(), center=Vector3()):
        self.material = material
        self.center = center

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
        self._radius = val if val > 0 else 0
    
    def get_json(self):
        return {'type' : 'sphere',
                'center' : self.center.get_json(),
                'radius' : self.radius,
                'medium' : self.material.get_json()}
    
class Block(GeometricObject):

    def __init__(self, size, **kwargs):
        self.size = size
        super().__init__(**kwargs)
    
    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, val):
        self._size = val if (val.x > 0 and val.y > 0 and val.z > 0) else Vector3()

    def get_json(self):
        return {"type" : "box", 
            "center" : self.center.get_json(), 
            "size" : self.size.get_json(),
            "medium" : self.material.get_json()}
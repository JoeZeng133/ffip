import math
import numbers
from math import cos, sin
from numbers import Number
from copy import deepcopy
import numpy as np
from abc import ABCMeta, abstractmethod

class Vector3(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = float(x) if type(x) is int else x
        self.y = float(y) if type(y) is int else y
        self.z = float(z) if type(z) is int else z

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z

        return Vector3(x, y, z)

    def __sub__(self, other):
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

class Medium:

    def __init__(self,
                 epsilon=1,
                 mu=1,
                 E_conductivity=0,
                 H_conductivity=0,
                 E_susceptibilities=[],
                 H_susceptibilities=[]):

        self.epsilon = epsilon
        self.mu = mu
        self.E_susceptibilities = E_susceptibilities
        self.H_susceptibilities = H_susceptibilities
        self.E_conductivity = E_conductivity
        self.H_conductivity = H_conductivity
        self.E_susceptibilities_array = None
        self.H_susceptibilities_array = None

    def get_suceptibilities_array(self, E_susceptibilities=[], H_susceptibilities=[]):
        lenE = len(E_susceptibilities)
        lenH = len(H_susceptibilities)

        if (self.E_susceptibilities_array is None and lenE > 0):
            self.E_susceptibilities_array = np.zeros(lenE)
            for esus in self.E_susceptibilities:
                for (i, item) in enumerate(E_susceptibilities):
                    if (esus == item):
                        self.E_susceptibilities_array[i] = esus.sigma
                        break


        if (self.H_susceptibilities_array is None and lenH > 0):
            self.H_susceptibilities_array = np.zeros(lenH)
            for hsus in self.H_susceptibilities:
                for (i, item) in enumerate(H_susceptibilities):
                    if (hsus == item):
                        self.E_susceptibilities_array[i] = esus.sigma
                        break

        return (self.E_susceptibilities_array, self.H_susceptibilities_array)
            

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
                

class LorentzianSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super.__init__(**kwargs)
        super.build_eqs(a0=frequency*frequency, b0=frequency*frequency, b1=gamma, b2=1)

class DrudeSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super.__init__(**kwargs)
        super.build_eqs(a0=frequency*frequency, b1=gamma, b2=1)

class CriticalPointsSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, phi=0.0, **kwargs):
        super.__init__(**kwargs)
        super.build_eqs(a0=2*frequency(frequency*cos(phi) - gamma*sin(phi)),
                        a1=-2*frequency*sin(phi),
                        b0=frequency*frequency+gamma*gamma,
                        b1=2*gamma,
                        b2=1)

class GeometricObject(object):

    def __init__(self, material=Medium(), center=Vector3(), esus=None, hsus=None):
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
        super.__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, val):
        self._radius = val if val > 0 else 0
    
class Block(GeometricObject):

    def __init__(self, size, **kwargs):
        self.size = size
        super.__init__(**kwargs)
    
    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, val):
        self._size = val if (val.x > 0 and val.y > 0 and val.z > 0) else Vector3()
from ffip.geom import Vector3, cmp_shape
from numpy import pi, exp
import numpy as np
from copy import deepcopy

class Source_Function:
    # functions for use in defining excitations
    def __init__(self):
        pass

    @property
    def start_time(self):
        return self._start_time
    
    @start_time.setter
    def start_time(self, val):
        self._start_time = val
    
    @property
    def end_time(self):
        return self._end_time
    
    @end_time.setter
    def end_time(self, val):
        self._end_time = val
    
    @property
    def duration(self):
        return self.end_time - self.start_time

class Gaussian1(Source_Function):
    # Gaussian first derivative

    def __init__(self, frequency, cutoff=4, start_time=0):
        
        self.frequency = float(frequency)
        self.cutoff = float(cutoff)
        self.start_time = float(start_time)
        
        # derived parameters
        self.width = 1 / (2 * pi * frequency)
        self.end_time = self.start_time + self.width * cutoff * 2
        self.offset_time = (self.start_time + self.end_time) / 2
    
    def get_json(self):
        return {'type' : 'Gaussian1',
                'frequency' : self.frequency,
                'cutoff' : self.cutoff,
                'start time' : self.start_time}
    
    def __call__(self, t):
        arg = (t - self.offset_time) / self.width
        return arg * exp(-(arg**2) / 2)

class Source:
    def __init__(self, function, center=Vector3(), size=Vector3(), amplitude=1, field_component='Ex'):
        self.function = deepcopy(function)
        self.center = center.copy()
        self.size = size.copy()
        self.amplitude = float(amplitude)
        self.field_component = deepcopy(field_component)

    def dump(self, file):
        pass

    def get_json(self):
        return {'type' : 'volume source',
                'center' : self.center.get_json(),
                'size' : self.size.get_json(),
                'amplitude' : self.amplitude,
                'field component' : self.field_component,
                'function' : self.function.get_json()}

class Inhom_Source:
    def __init__(self, function, amplitude=None, frequency=1, center=Vector3(), size=Vector3(), dim=Vector3(), field_component='Ex', suffix='0'):
        self.dim = dim.round()
        self.function = deepcopy(function)
        self.frequency = float(frequency)
        self.center = center.copy()
        self.size = size.copy()
        self.field_component = deepcopy(field_component)

        if callable(amplitude):
            self.amplitude = amplitude(np.stack(np.meshgrid(self.z, self.y, self.x, indexing='ij'), axis=-1))

        elif isinstance(amplitude, np.ndarray):
            self.amplitude = amplitude

        elif amplitude is None:
            self.amplitude = np.zeros(self.shape, dtype=complex)
        else:
            raise ValueError('amplitude input is not numpy array or a function')
        
        self.amplitude_dataset_name = 'inhom source amplitude %s' % suffix
        self.phase_dataset_name = 'inhom source phase %s' % suffix
    
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
    def amplitude(self):
        return self._amp   
    
    @property
    def shape(self):
        return (int(self.dim.z), int(self.dim.y), int(self.dim.x))
    
    @amplitude.setter
    def amplitude(self, val):
        if not cmp_shape(self.shape, val.shape):
            raise ValueError('amplitude shape is not compatibale with source dimension')
        self._amp = np.array(val, copy=True)
    
    @property
    def numel(self):
        return int(self.dim.prod())

    def dump(self, file):
        file.create_dataset(self.amplitude_dataset_name, data=np.abs(self.amplitude.ravel()))
        file.create_dataset(self.phase_dataset_name, data=np.angle(self.amplitude.ravel()))
    
    def get_json(self):
        return {'type' : 'inhomogeneous source',
                'center' : self.center.get_json(),
                'size' : self.size.get_json(),
                'dimension' : self.dim.get_json(),
                'field component' : self.field_component,
                'amplitude dataset' : self.amplitude_dataset_name,
                'phase dataset' : self.phase_dataset_name,
                'frequency' : self.frequency,
                'function' : self.function.get_json()}

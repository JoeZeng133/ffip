from ffip.geom import Vector3
from numpy import pi, exp

class Source_Function:
    # functions for use in defining excitations
    def __init__(self, start_time=0, end_time=1e20):
        self.start_time = start_time
        self.end_time = end_time

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


class Gaussian1(Source_Function):
    # Gaussian first derivative

    def __init__(self, frequency, cutoff=4, start_time=0):
        
        self.frequency = frequency
        self.cutoff = cutoff
        self.start_time = start_time
        
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
        self.function = function
        self.center = center
        self.size = size
        self.amplitude = float(amplitude)
        self.field_component = field_component

    def get_json(self):
        return {'type' : 'volume source',
                'center' : self.center.get_json(),
                'size' : self.size.get_json(),
                'amplitude' : self.amplitude,
                'field component' : self.field_component,
                'function' : self.function.get_json()}
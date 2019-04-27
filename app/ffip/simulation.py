import ffip
import numpy as np
import json
import h5py
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import trapz

class Fields_DFT:

    def __init__(self, center, size, frequency, field_component, count=0):
        self.count = count
        self.group_name = 'volume fields dft' + str(self.count)
        self.center = center
        self.size = size
        
        self.frequency = np.sort(np.asarray(frequency, dtype=float))
        self.field_component = field_component

    def get_json(self):
        return {'type' : 'volume fields dft',
                'center' : self.center.get_json(),
                'size' : self.size.get_json(),
                'frequency' : self.frequency.tolist(),
                'field component' : self.field_component,
                'output group' : self.group_name}

    @property
    def values(self):
        return self.v

    def read(self, h5file):
        group = h5file[self.group_name]
        self.x = group['x'][:]
        self.y = group['y'][:]
        self.z = group['z'][:]
        self.v = group['real'][:] + 1j * group['imag'][:]
    
    def get_interpolant(self, substract=None, **kwargs):
        if substract is not None and not isinstance(substract, Fields_DFT):
            raise TypeError('Substract is not an instance of Fields_DFT')
            
        return RegularGridInterpolator(
            (self.frequency, self.z, self.y, self.x), 
            self.values - substract.values if substract is not None else self.values, 
            **kwargs)
        
    def get_metadata(self):
        return self.z, self.y, self.x

class Flux_Region:
    def __init__(self, sim, center, size, frequency, weight):
        dim = (size.x != 0) + (size.y != 0) + (size.z != 0)
        if (dim != 2):
            raise NameError('Flux region dimension %d is invalid' % dim)
        
        self.direction = 0 if size.x == 0 else 1 if size.y == 0 else 2
        self.weight = weight
        self.center = center
        self.size = size
        self.frequency = np.sort(np.asarray(frequency))
        self.dx = sim.dx

        dict = {0 : 'x', 1 : 'y', 2 : 'z'}
        e1 = 'E' + dict[(self.direction + 1) % 3]
        h2 = 'H' + dict[(self.direction + 2) % 3]
        self._e1_dft = sim.add_dft_fields(center=center, size=size, frequency=frequency, field_component=e1)
        self._h2_dft = sim.add_dft_fields(center=center, size=size, frequency=frequency, field_component=h2)

    @property
    def e1_dft(self):
        return self._e1_dft
    
    @property
    def h2_dft(self):
        return self._h2_dft

    def get_values(self, scale=1, substract=None):
        if substract is not None and not isinstance(substract, Flux_Region):
            raise TypeError('Substract is not an instance of Flux_Region')
 
        len = (self.size / self.dx * scale).round() + ffip.Vector3(1,1,1)
        p1 = self.center - self.size / 2
        p2 = self.center + self.size / 2

        x, dx = np.linspace(p1.x, p2.x, len.x, retstep=True)
        y, dy = np.linspace(p1.y, p2.y, len.y, retstep=True)
        z, dz = np.linspace(p1.z, p2.z, len.z, retstep=True)
        pts = np.stack(np.meshgrid(self.frequency, z, y, x, indexing='ij'), axis=-1)

        # consider degenerating geometry
        meas = 1
        if not np.isnan(dx):
            meas *= dx
        if not np.isnan(dy):
            meas *= dy
        if not np.isnan(dz):
            meas *= dz

        e1_interp = self.e1_dft.get_interpolant(
            substract=substract.e1_dft if substract is not None else None, 
            method='nearest',bounds_error=False, fill_value=None)
        
        h2_interp = self.h2_dft.get_interpolant(
            substract=substract.h2_dft if substract is not None else None, 
            method='nearest',bounds_error=False, fill_value=None)

        e1 = e1_interp(pts)
        h2 = h2_interp(pts)
        flux = np.real(e1 * np.conj(h2)) * meas

        # integration along dimensions (1, 2, 3) = (dz, dy, dx)
        for i in range(3):
            flux = trapz(flux, axis=1) if flux.shape[1] > 1 else np.squeeze(flux, axis=1)

        return flux * self.weight

    def get_interpolant(self, scale=1, substract=None):
        flux = self.get_values(scale=scale, substract=substract)
        return RegularGridInterpolator((self.frequency,), flux)

class Flux_Box:
    def __init__(self, sim, center, size, frequency):
        self._flux_regions = []
        self.frequency = np.sort(np.asarray(frequency))

        for dir in range(3):
            for side in [-1, 1]:
                p1 = center - size/2
                p2 = center + size/2
                if side == -1:
                    p2[dir] = p1[dir]
                else:
                    p1[dir] = p2[dir]
                self._flux_regions.append(sim.add_flux_region(center=center, size=size, frequency=frequency, weight=side))
    
    @property
    def flux_regions(self):
        return self._flux_regions

    def get_values(self, scale=1, substract=None):
        if substract is not None and not isinstance(substract, Flux_Box):
            raise TypeError('Substract is not an instance of Flux_Box')
        
        flux = 0
        for i in range(6):
            flux += self.flux_regions[i].get_values(
                scale=scale, 
                substract=substract.flux_regions[i] if substract is not None else None)
        
        return flux
        
    def get_interpolant(self, scale=1, substract=None):
        flux = self.get_values(scale=scale, substract=substract)
        return RegularGridInterpolator((self.frequency), flux)

class Simulation:
    def __init__(self,
                 size,
                 resolution,
                 pmls = [],
                 geometry=[],
                 sources=[],
                 symmetry=[],
                 default_material=ffip.Medium(),
                 periodic=False,
                 progress_interval=4,
                 courant=0.5,
                 input_file='empty.h5',
                 fields_output_file='output.h5'):
                 
            print("Simulation created\n")
            self.size = size
            self.resolution = resolution
            self.geometry = geometry
            self.sources = sources
            self.symmetry = symmetry
            self.default_material = default_material
            self.periodic = periodic
            self.progress_interval = progress_interval
            self.courant = courant
            self.pmls = pmls
            self.input_file = input_file
            self.fields_output_file = fields_output_file

            # post processing
            self.dft_fields = []
            self.flux_regions = []
            self.flux_box = []

            # derive sigma_max_opt
            for pml in self.pmls:
                if pml.sigma_max is None:
                    pml.sigma_max = 0.8 * (1 + pml.order) * resolution
    @property
    def dx(self):
        return 1/self.resolution

    def get_json(self, stop_condition):
        res = {}
        # required parameters
        res['courant number'] = self.courant
        res['resolution'] = self.resolution
        res['size'] = self.size.get_json()
        res['default material'] = self.default_material.get_json()
        res['progress interval'] = self.progress_interval
        res['input file'] = self.input_file
        res['fields output file'] = self.fields_output_file
        res['stop condition'] = stop_condition.get_json()
        
        # optional parameters
        if self.symmetry:
            res['symmetry'] = [item.get_json() for item in self.periodic]
        
        if self.pmls:
            res['PML'] = [item.get_json() for item in self.pmls]
        
        if self.periodic:
            res['periodic'] = self.periodic
        
        if self.geometry:
            res['geometry'] = [item.get_json() for item in self.geometry]
        
        if self.dft_fields:
            res['fields output'] = [item.get_json() for item in self.dft_fields]
        
        if self.sources:
            res['sources'] = [item.get_json() for item in self.sources]
        
        return res
    
    def run(self, stop_condition):
        json.dump(self.get_json(stop_condition), open('config.json', 'w'), indent=4)
        input('Press Enter when Simulation finished')
        output_file = h5py.File(self.fields_output_file, 'r')
        for item in self.dft_fields:
            item.read(output_file)
    
    def add_dft_fields(self, center=ffip.Vector3(), size=ffip.Vector3(), frequency=[], field_component='Ex'):
        res = Fields_DFT(center, size, frequency, field_component, len(self.dft_fields))
        self.dft_fields.append(res)
        return res

    def add_flux_region(self, center=ffip.Vector3(), size=ffip.Vector3(), frequency=[], weight=1):
        res = Flux_Region(self, center=center, size=size, frequency=frequency, weight=weight)
        self.flux_regions.append(res)
        return res
    
    def add_flux_box(self, center=ffip.Vector3(), size=ffip.Vector3(), frequency=[]):
        res = Flux_Box(self, center=center, size=size, frequency=frequency)
        self.flux_box.append(res)
        return res

class Symmetry:
    def __init__(self, direction='x', phase_factor=1):
        self.direction = direction
        self.phase_factor = int(phase_factor)
    
    def get_json(self):
        return {'direction' : self.direction,
                'phase factor' : self.phase_factor}

class run_until_time:
    def __init__(self, time=0):
        self.time = float(time)
    
    def get_json(self):
        return {'type' : 'time', 
                'time' : self.time}

class run_until_fields_decay:
    def __init__(self, position=ffip.Vector3(), field_component='Ez', time_interval_examined=1, decayed_by=1e-3):
        self.position = position
        self.field_component = field_component
        self.time_interval_examined = float(time_interval_examined)
        self.decayed_by = float(decayed_by)

    def get_json(self):
        return {'type' : 'decay',
                'position' : self.position.get_json(),
                'field component' : self.field_component,
                'time interval examined' : self.time_interval_examined,
                'decayed by' : self.decayed_by}

class PML:
    def __init__(self, thickness, direction=None, side=None, k_max=1, order=3, sigma_max=None):
        self.thickness = float(thickness)
        self.direction = direction
        self.side = side
        self.k_max = float(k_max)
        self.order = int(order)
        self.sigma_max = None
    

    def get_json(self):
        res = {}
        res['thickness'] = self.thickness

        if self.direction is not None:
            res['direction'] = self.direction
        
        if self.side is not None:
            res['side'] = self.side
        
        res['k max'] = self.k_max if self.k_max is not None else 1.0
        res['order'] = self.order if self.order is not None else 3
        res['sigma max'] = self.sigma_max

        return res

class Flux_region:
    pass
    




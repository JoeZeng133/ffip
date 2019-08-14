import subprocess
import ffip
import numpy as np
import json
import h5py
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import trapz
from ffip.geom import Medium, Vector3, Two_Medium_Box, General_Medium_Box,cmp_shape, metadata_to_pts
from math import pi
from ffip.source import Source, Inhom_Source
from copy import deepcopy, copy


direction_dict = {'x': 0, 'y': 1, 'z': 2}
side_dict = {'positive': 1, 'negative': -1}

class Fields_DFT:

    def __init__(self, center, size, frequency, field_component, suffix='0', degenerate=False, dx=0):
        self.suffix = str(suffix)
        self.group_name = 'volume fields dft %s' % suffix
        self.center = center.copy()
        self.size = size.copy()
        self.degenerate = bool(degenerate)
        self.dx = float(dx)

        if not self.degenerate and dx == 0:
            raise ValueError('for non degenerate dft, dx must be provided')

        self.frequency = np.sort(np.asarray(frequency, dtype=float))
        # make a copy of original frequency
        self.frequency_config = self.frequency.tolist()
        self.field_component = field_component

    def get_json(self):
        return {'type': 'volume fields dft',
                'center': self.center.get_json(),
                'size': self.size.get_json(),
                'frequency': self.frequency_config,
                'field component': self.field_component,
                'output group': self.group_name}

    @property
    def values(self):
        return self.v

    def read(self, h5file):
        group = h5file[self.group_name]
        self.frequency = group['f'][:]
        self.x = group['x'][:]
        self.y = group['y'][:]
        self.z = group['z'][:]
        self.v = group['real'][:] + 1j * group['imag'][:]

        if not self.degenerate:
            # this is to cope with interpolation not working with dimension size 1
            # expanding frequency is different from expanding coordinates
            if self.frequency.size == 1:
                self.frequency = self.frequency * np.array([0.9, 1.1])
            
            if self.x.size == 1:
                self.x = self.x + np.array([-self.dx/2, self.dx/2])

            if self.y.size == 1:
                self.y = self.y + np.array([-self.dx/2, self.dx/2])
            
            if self.z.size == 1:
                self.z = self.z + np.array([-self.dx/2, self.dx/2])

            for i in range(4):
                if self.v.shape[i] == 1:
                    self.v = np.concatenate((self.v, self.v), axis=i)
            
    
    def __call__(self, f, z, y, x, substract=None, **kwargs):

        interp = self.get_interpolant(substract, **kwargs)
        pts = np.stack(np.meshgrid(f, z, y, x, indexing='ij'), axis=-1)
        return interp(pts)

    def get_interpolant(self, substract=None, **kwargs):
        if substract is not None and not isinstance(substract, Fields_DFT):
            raise TypeError('Substract is not an instance of Fields_DFT')

        return RegularGridInterpolator(
            (self.frequency, self.z, self.y, self.x),
            self.values - substract.values if substract is not None else self.values,
            **kwargs)

    def get_metadata(self):
        return self.frequency, self.z, self.y, self.x


class Flux_Region:
    def __init__(self, sim, center, size, frequency, weight):
        dim = (size.x != 0) + (size.y != 0) + (size.z != 0)
        if (dim != 2):
            raise NameError('Flux region dimension %d is invalid' % dim)

        self.direction = 0 if size.x == 0 else 1 if size.y == 0 else 2
        self.weight = np.array(weight, copy=True)
        self.center = center.copy()
        self.size = size.copy()
        self.frequency = np.sort(np.asarray(frequency))
        self.dx = sim.dx

        self._e1_dft = sim.add_dft_fields(
            center=center, size=size, frequency=frequency, field_component=self.e1_str)
        self._h2_dft = sim.add_dft_fields(
            center=center, size=size, frequency=frequency, field_component=self.h2_str)
    
    @property
    def e1_str(self):
        dict = {0: 'x', 1: 'y', 2: 'z'}
        return 'E' + dict[(self.direction + 1) % 3]
    
    @property
    def h2_str(self):
        dict = {0: 'x', 1: 'y', 2: 'z'}
        return 'H' + dict[(self.direction + 2) % 3]

    @property
    def e1_dft(self):
        return self._e1_dft

    @property
    def h2_dft(self):
        return self._h2_dft

    def get_fields(self, scale=1, substract=None, **kwargs):
        # return e1, h2
        if substract is not None and not isinstance(substract, Flux_Region):
            raise TypeError('Substract is not an instance of Flux_Region')

        e1_interp = self.e1_dft.get_interpolant(
            substract=substract.e1_dft if substract is not None else None,
            **kwargs)

        h2_interp = self.h2_dft.get_interpolant(
            substract=substract.h2_dft if substract is not None else None,
            **kwargs)
        
        pts = self.get_pts(scale=scale)
        e1 = e1_interp(pts)
        h2 = h2_interp(pts)

        return e1, h2
    
    def get_len(self, scale=1):
        return (self.size / self.dx * scale).round() + ffip.Vector3(1, 1, 1)
    
    def get_metadata(self, scale=1):
        len = self.get_len(scale=scale)
        p1 = self.center - self.size / 2
        p2 = self.center + self.size / 2

        x = np.linspace(p1.x, p2.x, len.x)
        y = np.linspace(p1.y, p2.y, len.y)
        z = np.linspace(p1.z, p2.z, len.z)

        return self.frequency, z, y, x

    def get_pts(self, scale=1):
        # return integration points

        return np.stack(np.meshgrid(*self.get_metadata(scale=scale), indexing='ij'), axis=-1)
    
    def get_int_weights(self, scale=1):
        # return integration weights

        len = self.get_len(scale=scale)

        x = np.ones((int(len.x),))
        y = np.ones((int(len.y),))
        z = np.ones((int(len.z),))

        meas = 1
        if len.x > 1:
            meas *= self.size.x / (len.x - 1)
        if len.y > 1:
            meas *= self.size.y / (len.y - 1)
        if len.z > 1:
            meas *= self.size.z / (len.z - 1)

        for item in [x, y, z]:
            if item.size > 1:
                item[0] = 0.5
                item[-1] = 0.5
        
        Z, Y, X = np.meshgrid(z, y, x, indexing='ij')
        return X * Y * Z * meas * self.weight
        
    def get_values(self, scale=1, substract=None, **kwargs):
        # return fluxes

        e1, h2 = self.get_fields(scale=scale, substract=substract, **kwargs)
        int_weights = self.get_int_weights(scale=scale)
        
        flux = np.real(e1 * np.conj(h2)) * int_weights

        # integration along dimensions (1, 2, 3) = (dz, dy, dx)
        flux = np.sum(flux, axis=(1, 2, 3))

        return flux

    def get_interpolant(self, scale=1, substract=None):
        # return fluxes as an interpolant

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
                self._flux_regions.append(sim.add_flux_region(
                    center=(p1+p2)/2, size=p2-p1, frequency=frequency, weight=side))

    @property
    def flux_regions(self):
        return self._flux_regions

    def get_flux_region(self, direction='x', side='negative'):
        return self.flux_regions[direction_dict[direction] * 2 + (0 if side_dict[side] == -1 else 1)]

    def get_values(self, scale=1, substract=None, **kwargs):
        if substract is not None and not isinstance(substract, Flux_Box):
            raise TypeError('Substract is not an instance of Flux_Box')

        flux = 0
        for i in range(6):
            flux += self.flux_regions[i].get_values(
                scale=scale,
                substract=substract.flux_regions[i] if substract is not None else None,
                **kwargs)

        return flux

    def get_interpolant(self, scale=1, substract=None, **kwargs):
        flux = self.get_values(scale=scale, substract=substract, **kwargs)
        return RegularGridInterpolator((self.frequency), flux)


class Geom_Info:
    def __init__(self, index=0, field_component='Ex'):
        self.field_component = str(field_component)
        self.index = int(index)
        self.gname = 'geometry info request %d %s' % (
            self.index, self.field_component)

    def read(self, file):
        group = file[self.gname]
        self._x = group['x'][:]
        self._y = group['y'][:]
        self._z = group['z'][:]

    def get_json(self):
        return {
            'field component': self.field_component,
            'geometry index': self.index,
            'output group': self.gname
        }

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z


class Simulation:
    def __init__(self,
                 size,
                 resolution,
                 pmls=[],
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
        self.size = size.copy()
        self.resolution = float(resolution)
        self.geometry = deepcopy(geometry)
        self.sources = deepcopy(sources)
        self.symmetry = deepcopy(symmetry)
        self.default_material = deepcopy(default_material)
        self.periodic = deepcopy(periodic)
        self.progress_interval = float(progress_interval)
        self.courant = float(courant)
        self.pmls = deepcopy(pmls)
        self.input_file = copy(input_file)
        self.fields_output_file = copy(fields_output_file)

        # post processing
        self.dft_fields = []
        self.flux_regions = []
        self.flux_boxes = []
        self.geom_info_requests = []

        # derive sigma_max_opt
        for pml in self.pmls:
            if pml.sigma_max is None:
                pml.sigma_max = 0.8 * (1 + pml.order) * resolution
    
    @property
    def input_file(self):
        return self._input_file
    
    @input_file.setter
    def input_file(self, val):
        if not isinstance(val, str):
            raise ValueError('file name is not a string')
        
        self._input_file = val
    
    @property
    def fields_output_file(self):
        return self._fields_output_file
    
    @fields_output_file.setter
    def fields_output_file(self, val):
        if not isinstance(val, str):
            raise ValueError('file name is not a string')

        self._fields_output_file = val

    @property
    def dx(self):
        return 1/self.resolution

    def dump(self, file):
        # dump input data into input_file
        if self.geometry:
            for geom in self.geometry:
                geom.dump(file)

        if self.sources:
            for src in self.sources:
                src.dump(file)

    def get_json(self, stop_condition, decomposition=None):
        # get configuration json strings
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
        if decomposition is not None:
            res['decomposition'] = decomposition

        if self.symmetry:
            res['symmetry'] = [item.get_json() for item in self.periodic]

        if self.pmls:
            res['PML'] = [item.get_json() for item in self.pmls]

        if self.periodic:
            res['periodic'] = self.periodic

        if self.geometry:
            res['geometry'] = [item.get_json() for item in self.geometry]

        if self.sources:
            res['sources'] = [item.get_json() for item in self.sources]

        if self.dft_fields:
            res['fields output'] = [item.get_json()
                                    for item in self.dft_fields]

        if self.geom_info_requests:
            res['geometry info request'] = [item.get_json()
                                            for item in self.geom_info_requests]

        return res

    def run(self, stop_condition, decomposition=None, np=1, skip=False, pop=False, **kwargs):

        if decomposition is not None:
            np = int(decomposition[0] * decomposition[1] * decomposition[2])

        with h5py.File(self.input_file, 'w') as input_file_h5:
            # dump json configuration file
            json.dump(self.get_json(stop_condition=stop_condition, decomposition=decomposition),
                    open('config.json', 'w'), indent=4)
            # dump input data
            self.dump(input_file_h5)

        # invoking externel bash command
        if not skip:
            process = subprocess.Popen('mpirun -np %d run_sim_json' % np, shell=True, **kwargs)
            process.wait()

        # run externel command manually
        if pop:
            input('press enter when finished')
        
        with h5py.File(self.fields_output_file, 'r') as output_file_h5:
            # read dft fields
            for item in self.dft_fields:
                item.read(output_file_h5)
            # read geometry info requests
            for item in self.geom_info_requests:
                item.read(output_file_h5)

    def add_geometry(self, geom):
        self.geometry.insert(0, geom)

    def add_source(self, src):
        self.sources.insert(0, src)

    def request_geom_info(self, geom, field_component='Ex'):
        # make a request for coordinates of (ew, hw) points inside a geometry
        for i in range(len(self.geometry)):
            if self.geometry[i] is geom:
                req = Geom_Info(index=i, field_component=field_component)
                self.geom_info_requests.append(req)
                return req

        raise ValueError('Geometry non-existed')

    def add_dft_fields(self, center=ffip.Vector3(), size=ffip.Vector3(), frequency=[], field_component='Ex'):
        res = Fields_DFT(center, size, frequency,
                         field_component, str(len(self.dft_fields)), degenerate=False, dx=self.dx)
        self.dft_fields.append(res)
        return res

    def add_flux_region(self, center=ffip.Vector3(), size=ffip.Vector3(), frequency=[], weight=1):
        res = Flux_Region(self, center=center, size=size,
                          frequency=frequency, weight=weight)
        self.flux_regions.append(res)
        return res

    def add_flux_box(self, center=ffip.Vector3(), size=ffip.Vector3(), frequency=[]):
        res = Flux_Box(self, center=center, size=size, frequency=frequency)
        self.flux_boxes.append(res)
        return res


class Symmetry:
    def __init__(self, direction='x', phase_factor=1):
        self.direction = direction
        self.phase_factor = int(phase_factor)

    def get_json(self):
        return {'direction': self.direction,
                'phase factor': self.phase_factor}


class run_until_time:
    def __init__(self, time=0):
        self.time = float(time)

    def get_json(self):
        return {'type': 'time',
                'time': self.time}
class run_until_dft:
    def __init__(self, center=ffip.Vector3(), size=ffip.Vector3(), frequency=1, field_component='Ex', time_interval_examined=1, var=1e-2):
        self.center = center.copy()
        self.size = size.copy()
        self.field_component = str(field_component)
        self.time_interval_examined = float(time_interval_examined)
        self.var = float(var)
        self.frequency = float(frequency)
    
    def get_json(self):
        return {'type' : 'dft',
                'center' : self.center.get_json(),
                'size'   : self.size.get_json(),
                'field component' : self.field_component,
                'time interval examined': self.time_interval_examined,
                'variation' : self.var,
                'frequency' : self.frequency
                }


class run_until_fields_decay:
    def __init__(self, position=ffip.Vector3(), field_component='Ez', time_interval_examined=1, decayed_by=1e-3):
        self.position = position.copy()
        self.field_component = str(field_component)
        self.time_interval_examined = float(time_interval_examined)
        self.decayed_by = float(decayed_by)

    def get_json(self):
        return {'type': 'decay',
                'position': self.position.get_json(),
                'field component': self.field_component,
                'time interval examined': self.time_interval_examined,
                'decayed by': self.decayed_by}


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

def field_match_objective(obj_val, component='Ex', weight=1):
    # return triplets (list of components, list of fprimes, objective function)
    # for adjoint source calculations

    if component in ['Ex', 'Ey', 'Ez']:
        # original complex component l2 norm
        def fun1(val):
            return np.sum(0.5 * np.ravel(weight) * np.abs(val.ravel() - obj_val.ravel())**2)
        
        def fpconj1(val):
            return (np.conj(val - obj_val) * weight,)
        
        return [component], fpconj1, fun1
        
    if component in ['|Ex|', '|Ey|', '|Ez|']:
        # absolute value l2 norm
        def fun2(val):
            return np.sum(0.5 * np.ravel(weight) * (np.abs(val).ravel() - obj_val.ravel())**2)
        
        def fpconj2(val):
            return ((1 - obj_val / np.abs(val)) * np.conj(val) * weight,)
        
        return [component[1:-1]], fpconj2, fun2
    
    if component == '|E|':
        # absolute value of vector l2 norm

        def fun3(ex, ey, ez):
            val = np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2)
            return np.sum(0.5 * (np.ravel(val) - np.ravel(obj_val))**2 * np.ravel(weight))
        
        def fpconj3(ex, ey, ez):
            val = np.sqrt(np.abs(ex)**2 + np.abs(ey)**2 + np.abs(ez)**2)
            return ((1 - obj_val / val) * np.conj(ex) * weight, 
                    (1 - obj_val / val) * np.conj(ey) * weight,
                    (1 - obj_val / val) * np.conj(ez) * weight
            )

        return ['Ex', 'Ey', 'Ez'], fpconj3, fun3


class Adjoint_Flux:
    def __init__(self, adjoint_simulation, forward_simulation, function, fluxes=[], frequency=1.0, scale=1):
        self.adjoint_simulation = adjoint_simulation
        self.forward_simulation = forward_simulation
        self.frequency = float(frequency)
        self.scale = float(scale)
        self.fluxes = fluxes
        self.e_adjoint_sources = []
        self.m_adjoint_sources = []


        for i in range(len(fluxes)):    
            flux = fluxes[i]

            if flux.frequency.size > 1:
                raise ValueError('Fluxes used in Adjoint_Flux cannot have more than 1 frequency')

            # j1 source
            self.e_adjoint_sources.append(
                Inhom_Source(
                        function=function,
                        frequency=frequency,
                        amplitude=None,
                        center=flux.center,
                        size=flux.size,
                        dim=flux.get_len(scale=scale),
                        field_component=flux.e1_str,
                        suffix='flux %d %s adjoint' % (i, flux.e1_str)
                    )
                )
            
            # m2 source
            self.m_adjoint_sources.append(
                Inhom_Source(
                        function=function,
                        frequency=frequency,
                        amplitude=None,
                        center=flux.center,
                        size=flux.size,
                        dim=flux.get_len(scale=scale),
                        field_component=flux.h2_str,
                        suffix='flux %d %s adjoint' % (i, flux.h2_str)
                    )
                )
            
            adjoint_simulation.add_source(self.e_adjoint_sources[-1])
            adjoint_simulation.add_source(self.m_adjoint_sources[-1])
    
    def eval_functionals_and_set_sources(self):
        res = 0

        for i in range(len(self.fluxes)):
            flux = self.fluxes[i]
            res += flux.get_values(scale=self.scale)
            e1, h2 = flux.get_fields(scale=self.scale)
            int_weights = flux.get_int_weights(scale=self.scale)

            self.e_adjoint_sources[i].amplitude = np.conj(h2[0, ...]) * int_weights
            self.m_adjoint_sources[i].amplitude = - np.conj(e1[0, ...]) * int_weights
        
        return res


class Adjoint_Source:
    comp_map = {'Ex': ['Ex'], 'Ey': ['Ey'], 'Ez': ['Ez'],
                '|Ex|': ['Ex'], '|Ey|': ['Ey'], '|Ez|': ['Ez'],
                '|E|': ['Ex', 'Ey', 'Ez']}

    def __init__(self, adjoint_simulation, forward_simulation, function, frequency=1.0, center=Vector3(), size=Vector3(), dim=Vector3(), functionals=[], udfs=[]):
        # list of functionals = [['ex', obj_val], ['|ey|', obj_val], ...]

        self.adjoint_simulation = adjoint_simulation
        self.forward_simulation = forward_simulation
        self.frequency = float(frequency)
        self.center = center.copy()
        self.size = size.copy()
        self.dim = dim.round()
        self.functionals = deepcopy(functionals)
        self.udfs = deepcopy(udfs)
        self.function = deepcopy(function)

        self.forward_dfts = {}
        self.adjoint_sources = {}
        self.forward_fields = {}

        p1 = self.center - self.size / 2
        self.x = np.linspace(p1.x, p1.x + self.size.x, dim.x)
        self.y = np.linspace(p1.y, p1.y + self.size.y, dim.y)
        self.z = np.linspace(p1.z, p1.z + self.size.z, dim.z)

        # legacy support
        for func in functionals:
            if not cmp_shape(self.shape, func[1].shape):
                raise ValueError('objective values does not match the shape')
            
            self.udfs.append(field_match_objective(func[1], func[0], weight=1))

        for udf in self.udfs:
            components = udf[0]
            fprime = udf[1]
            fun = udf[2]

            if not callable(fprime):
                raise ValueError('The derivative function is not callable')

            # if len(components) != len(fprimes):
            #     raise ValueError('The number of components does not match number of fprimes')
            
            if not callable(fun):
                raise ValueError('The objective function is not callable')

            for component in components:
                self.add_component(component)

    
    def add_component(self, component='Ex'):
        # add a field component into adjoint source and forward dfts
        
        if component not in self.forward_dfts:
            self.forward_dfts[component] = self.forward_simulation.add_dft_fields(
                center=self.center,
                size=self.size,
                frequency=[self.frequency],
                field_component=component
            )

            self.adjoint_sources[component] = Inhom_Source(
                function=self.function,
                frequency=self.frequency,
                amplitude=None,
                center=self.center,
                size=self.size,
                dim=self.dim,
                field_component=component,
                suffix=component + ' adjoint'
            )

            self.adjoint_simulation.add_source(self.adjoint_sources[component])

    @property
    def shape(self):
        return (int(self.dim.z), int(self.dim.y), int(self.dim.x))

    def eval_functionals_and_set_sources(self):

        # get forward fields
        for component in self.forward_dfts.keys():

            self.adjoint_sources[component].amplitude = np.zeros(self.shape, dtype=complex)
            self.forward_fields[component] = self.forward_dfts[component](
                    self.frequency, self.z, self.y, self.x).reshape(self.shape)

        res = 0
        for udf in self.udfs:
            components = udf[0]
            fprime = udf[1]
            fun = udf[2]

            args = []
            for component in components:
                args.append(self.forward_fields[component])
            
            res =  res + fun(*args)
            derivatives = fprime(*args)

            for i in range(len(components)):
                self.adjoint_sources[components[i]].amplitude += derivatives[i]
        
        return res

    @property
    def numel(self):
        return int(self.dim.prod())

class Adjoint_Volume:
    def __init__(self, adjoint_simulation, forward_simulation, frequency=1.0, center=Vector3(), size=Vector3(), dim=Vector3(), 
        density=None, medium1=Medium(), medium2=Medium(), 
        epsilon_fun=None, epsilon_der=None, e_sus = [], norm=1):
        # legacy support: e = e1 * density + (1 - density) * e2
        # general e = e(rho) = eps_inf + rho1 * e_sus1 + rho2 * e_sus2, de/drho provided

        self.adjoint_simulation = adjoint_simulation
        self.forward_simulation = forward_simulation
        self.frequency = float(frequency)
        self.center = center.copy()
        self.size = size.copy()
        self.dim = dim.round()
        self.norm = complex(norm)

        # use suffix adjoint to avoid dataset name conflicts
        if epsilon_fun is None:
            e1 = medium1.get_epsilon(self.frequency)
            e2 = medium2.get_epsilon(self.frequency)
            self.epsilon_fun = lambda rho : e1 * rho + (1 - rho) * e2
            self.epsilon_der = lambda rho : e1 - e2
            self.medium1 = deepcopy(medium1)
            self.medium2 = deepcopy(medium2)

            self.geom = Two_Medium_Box(
                size=size, 
                center=center, 
                dim=dim,
                density=density, 
                medium1=medium1, 
                medium2=medium2, 
                suffix='adjoint'
            )
            
        else:
            self.epsilon_der = deepcopy(epsilon_der)
            self.epsilon_fun = deepcopy(epsilon_fun)

            self.geom = General_Medium_Box(
                size=size,
                center=center,
                dim=dim,
                frequency=frequency,
                density=density,
                epsilon_fun=epsilon_fun,
                e_sus=e_sus,
                suffix='adjoint'
            )
        
        adjoint_simulation.add_geometry(self.geom)
        forward_simulation.add_geometry(self.geom)

        self.forward_dfts = [forward_simulation.add_dft_fields(
            center=center,
            size=size,
            frequency=[frequency],
            field_component=component
        ) for component in ['Ex', 'Ey', 'Ez']
        ]

        self.adjoint_dfts = [adjoint_simulation.add_dft_fields(
            center=center,
            size=size,
            frequency=[frequency],
            field_component=component
        ) for component in ['Ex', 'Ey', 'Ez']
        ]

        self.geom_infos = [forward_simulation.request_geom_info(
            geom=self.geom, 
            field_component=component
            ) for component in ['Ex', 'Ey', 'Ez']
        ]
    
    def forward_dft(self, field_component='Ex'):
        map = {'Ex' : 0, 'Ey' : 1, 'Ez' : 2}
        return self.forward_dfts[map[field_component]]
    
    def adjoint_dft(self, field_component='Ex'):
        map = {'Ex' : 0, 'Ey' : 1, 'Ez' : 2}
        return self.adjoint_dfts[map[field_component]]

    def get_sensitivity(self):
        self.pts = [np.stack((self.frequency * np.ones(item.z.shape), item.z, item.y, item.x), axis=-1) for item in self.geom_infos]

        rho_interp = RegularGridInterpolator(
            (self.z, self.y, self.x),
            self.density,
            bounds_error=True)
        
        self.rho = [rho_interp(np.stack((item.z, item.y, item.x), axis=-1)) for item in self.geom_infos]
        self.forward_fields = [self.forward_dfts[i].get_interpolant(method='linear', bounds_error=True)(self.pts[i]) for i in range(3)]
        self.adjoint_fields = [self.adjoint_dfts[i].get_interpolant(method='linear', bounds_error=True)(self.pts[i]) for i in range(3)]
        self.se = []
        self.se_transposed = []

        res = 0

        for i in range(3):
            self.se.append(np.real(2j * pi * self.frequency * self.forward_fields[i] * self.adjoint_fields[i] * self.epsilon_der(self.rho[i]) / self.norm))

            self.se_transposed.append(
                ffip.transpose(
                    self.x, self.y, self.z, 
                    np.stack((self.geom_infos[i].x, self.geom_infos[i].y, self.geom_infos[i].z,  self.se[i]), axis=-1)
                    )
            )

            res += self.se_transposed[i]

        return np.reshape(res, self.shape)
    
    def get_sensitivity2(self):
        pts = np.stack(np.meshgrid(self.frequency, self.z, self.y, self.x, indexing='ij'), axis=-1)
        forward_fields = [self.forward_dfts[i].get_interpolant(method='linear', bounds_error=True)(pts) for i in range(3)]
        adjoint_fields = [self.adjoint_dfts[i].get_interpolant(method='linear', bounds_error=True)(pts) for i in range(3)]
        res = 0

        for i in range(3):
            res = res + 2j * pi * self.frequency * forward_fields[i] * adjoint_fields[i] / self.norm
        
        res = np.reshape(res, self.shape)
        res = np.real(res * self.epsilon_der(self.density))
        
        return res

    @property
    def x(self):
        return self.geom.x
    
    @property
    def y(self):
        return self.geom.y
    
    @property
    def z(self):
        return self.geom.z

    @property
    def p1(self):
        return self.center - self.size/2

    @property
    def density(self):
        return self.geom.density
        
    @density.setter
    def density(self, val):
        self.geom.density = val
    
    @property
    def numel(self):
        return int(self.dim.prod())
    
    @property
    def shape(self):
        return self.geom.shape

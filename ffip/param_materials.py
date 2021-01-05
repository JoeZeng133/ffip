from ffip.geom import Param_Medium
from ffip.materials import *
import numpy as np

#------------------------------------------------------------------
# gold (Au) 500-1000 nm non-linear interpolation parameterized model
# based on IIR2 and IIR2 sqrt material
def Au_nonlin2():
	def e_inf_fun(rho: np.ndarray):
		return np.ones(np.shape(rho))

	def esus_amp_fun(rho: np.ndarray):
		sus1 = (rho**2 for i in range(len(Au_IIR2_susc)))
		sus2 = (2*rho*(1-rho) for i in range(len(Au_sqrt_IIR2_susc)))
		return np.stack((*sus1, *sus2), axis=-1)
	
	def e_inf_prime_fun(rho: np.ndarray):
		return np.zeros(np.shape(rho))
	
	def esus_amp_prime_fun(rho: np.ndarray):
		sus1 = (2 * rho for i in range(len(Au_IIR2_susc)))
		sus2 = ((2 - 4 * rho) for i in range(len(Au_sqrt_IIR2_susc)))
		return np.stack((*sus1, *sus2), axis=-1)

	m = Medium(epsilon=1, E_susceptibilities=Au_IIR2_susc + Au_sqrt_IIR2_susc)
	return Param_Medium(m, e_inf_fun, esus_amp_fun), Param_Medium(m, e_inf_prime_fun, esus_amp_prime_fun)

#------------------------------------------------------------------
# gold (Au) 500-1000 nm linear interpolation parameterized model
# based on IIR2 material
def Au_lin2():
	def e_inf_fun(rho: np.ndarray):
		return np.ones(np.shape(rho))

	def esus_amp_fun(rho: np.ndarray):
		sus1 = (rho for i in range(len(Au_IIR2_susc)))
		return np.stack((*sus1,), axis=-1)
	
	def e_inf_prime_fun(rho: np.ndarray):
		return np.zeros(np.shape(rho))
	
	def esus_amp_prime_fun(rho: np.ndarray):
		sus1 = (np.ones(np.shape(rho)) for i in range(len(Au_IIR2_susc)))
		return np.stack((*sus1,), axis=-1)

	m = Au_IIR2
	return Param_Medium(m, e_inf_fun, esus_amp_fun), Param_Medium(m, e_inf_prime_fun, esus_amp_prime_fun)

import ffip
import numpy as np

#%% test Au IIR2
m1 = ffip.Au
m2 = ffip.Au_IIR2
ft = 1/800

er1 = m1.get_epsilon(ft)
er2 = m2.get_epsilon(ft)

# print(er1)
# print(er2)

#%% test non-linear param medium
m,mp = ffip.Au_nonlin2()
rho = np.linspace(0, 1, 11)
ft = 1 / 800
er1 = m.get_dis_epsilon(rho, ft, 1)
er2 = (rho * np.sqrt(ffip.Au.get_dis_epsilon(ft, 1)) + (1-rho))**2

er1p = mp.get_dis_epsilon(rho, ft, 1)
er2p = 2 * (rho * np.sqrt(ffip.Au.get_dis_epsilon(ft, 1)) + (1-rho)) * (np.sqrt(ffip.Au.get_dis_epsilon(ft, 1)) - 1)

# print(er1)
# print(er2)

#%% test linear param medium
m,mp = ffip.Au_lin2()
rho = np.linspace(0, 1, 11)
ft = 1 / 800
er1 = m.get_dis_epsilon(rho, ft, 1)
er2 = ffip.Au.get_dis_epsilon(ft, 1) * rho + 1 - rho

er1p = mp.get_dis_epsilon(rho, ft, 1)
er2p = np.ones(np.shape(rho)) * (ffip.Au.get_dis_epsilon(ft, 1) - 1)

# print(m.esus_amp_fun(rho))
print(er1)
print(er2)
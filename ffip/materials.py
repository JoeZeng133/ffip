from ffip.geom import DrudeSusceptibility, LorentzianSusceptibility, Medium, DeybeSusceptibility, IIR_Susceptibility
from math import sqrt, pi
from numpy import real, imag
# import ffip

# default unit length is 1 um
um_scale = 1.0e-3

# conversion factor for eV to c/um [=1/hc]
eV_um_scale = um_scale/1.23984193
# conversion factor for hz to c/um
hz_scale = 1e-6 / 3e8 * um_scale


#------------------------------------------------------------------
# FePt (single frequency model using Lorentz Pole)
def FePt(frequency=1.0):

    e = (3.2 + 2.6j)**2

    eps2 = real(e)
    sus2 = LorentzianSusceptibility(
        frequency=frequency,
        gamma=frequency,
        sigma=imag(e)
    )

    return Medium(epsilon=eps2, E_susceptibilities=[sus2])

#------------------------------------------------------------------
# gold (Au) L4D model
Au_plasma_frq = 9.03*eV_um_scale
Au_f0 = 0.760
Au_frq0 = 1e-10*eV_um_scale*sqrt(6.1971084e+21)
Au_gam0 = 0.053*eV_um_scale
Au_sig0 = Au_f0*Au_plasma_frq**2/Au_frq0**2
Au_f1 = 0.024
Au_frq1 = 0.415*eV_um_scale      # 2.988 um
Au_gam1 = 0.241*eV_um_scale
Au_sig1 = Au_f1*Au_plasma_frq**2/Au_frq1**2
Au_f2 = 0.010
Au_frq2 = 0.830*eV_um_scale      # 1.494 um
Au_gam2 = 0.345*eV_um_scale
Au_sig2 = Au_f2*Au_plasma_frq**2/Au_frq2**2
Au_f3 = 0.071
Au_frq3 = 2.969*eV_um_scale      # 0.418 um
Au_gam3 = 0.870*eV_um_scale
Au_sig3 = Au_f3*Au_plasma_frq**2/Au_frq3**2
Au_f4 = 0.601
Au_frq4 = 4.304*eV_um_scale      # 0.288 um
Au_gam4 = 2.494*eV_um_scale
Au_sig4 = Au_f4*Au_plasma_frq**2/Au_frq4**2

Au_susc = [DrudeSusceptibility(frequency=Au_frq0, gamma=Au_gam0, sigma=Au_sig0),
           LorentzianSusceptibility(frequency=Au_frq1, gamma=Au_gam1, sigma=Au_sig1),
           LorentzianSusceptibility(frequency=Au_frq2, gamma=Au_gam2, sigma=Au_sig2),
           LorentzianSusceptibility(frequency=Au_frq3, gamma=Au_gam3, sigma=Au_sig3),
           LorentzianSusceptibility(frequency=Au_frq4, gamma=Au_gam4, sigma=Au_sig4)]

Au = Medium(epsilon=1.0, E_susceptibilities=Au_susc)

#------------------------------------------------------------------
# gold (Au)
# fit to P.B. Johnson and R.W. Christy, Physical Review B, Vol. 6, pp. 4370-9, 1972

Au_JC_visible_frq0 = 1/(0.139779231751333*um_scale)
Au_JC_visible_gam0 = 1/(26.1269913352870*um_scale)
Au_JC_visible_sig0 = 1

Au_JC_visible_frq1 = 1/(0.404064525036786*um_scale)
Au_JC_visible_gam1 = 1/(1.12834046202759*um_scale)
Au_JC_visible_sig1 = 2.07118534879440

Au_JC_visible_susc = [DrudeSusceptibility(frequency=Au_JC_visible_frq0, gamma=Au_JC_visible_gam0, sigma=Au_JC_visible_sig0),
                      LorentzianSusceptibility(frequency=Au_JC_visible_frq1, gamma=Au_JC_visible_gam1, sigma=Au_JC_visible_sig1)]

Au_JC_visible = Medium(epsilon=6.1599, E_susceptibilities=Au_JC_visible_susc)

#------------------------------------------------------------------
# gold (Au)
# fit to E.D. Palik, Handbook of Optical Constants, Academic Press, 1985 
Au_visible_frq0 = 1/(0.0473629248511456*um_scale)
Au_visible_gam0 = 1/(0.255476199605166*um_scale)
Au_visible_sig0 = 1

Au_visible_frq1 = 1/(0.800619321082804*um_scale)
Au_visible_gam1 = 1/(0.381870287531951*um_scale)
Au_visible_sig1 = -169.060953137985

Au_visible_susc = [DrudeSusceptibility(frequency=Au_visible_frq0, gamma=Au_visible_gam0, sigma=Au_visible_sig0),
                   LorentzianSusceptibility(frequency=Au_visible_frq1, gamma=Au_visible_gam1, sigma=Au_visible_sig1)]

Au_visible = Medium(epsilon=0.6888, E_susceptibilities=Au_visible_susc)

#------------------------------------------------------------------
# gold (Au) LD model
Au_LD_frq0 = 2113.6e12 * hz_scale
Au_LD_gam0 = 15.92e12 * hz_scale
Au_LD_sig0 = 1

Au_LD_frq1 = 650.07e12 * hz_scale
Au_LD_gam1 = 104.86e12 * hz_scale
Au_LD_sig1 = 1.09

Au_LD_susc = [DrudeSusceptibility(frequency=Au_LD_frq0, gamma=Au_LD_gam0, sigma=Au_LD_sig0),
            LorentzianSusceptibility(frequency=Au_LD_frq1, gamma=Au_LD_gam1, sigma=Au_LD_sig1)]
Au_LD = Medium(epsilon=5.9673, E_susceptibilities=Au_LD_susc)

#------------------------------------------------------------------
# gold (Au) IIR 200-1000 nm model
# Au_IIR1_susc = [
# 	IIR_Susceptibility((-1.1859e+15 + 5.1317e+15j) * hz_scale, (5.0833e+14 - 7.9338e+15j) * hz_scale),
# 	IIR_Susceptibility((-5.5450e+13 + 2.2427e+14j) * hz_scale, (-2.6683e+13 - 3.2983e+17j) * hz_scale)]
# Au_IIR1 = Medium(epsilon=1, E_susceptibilities=Au_IIR1_susc)

#------------------------------------------------------------------
# gold (Au) IIR 500-1000 nm model
Au_IIR2_susc = [
	IIR_Susceptibility((-1.1859e+15 + 5.1317e+15j) * hz_scale, (5.0833e+14 - 7.9338e+15j) * hz_scale),
	IIR_Susceptibility((-5.5450e+13 + 2.2427e+14j) * hz_scale, (-2.6683e+13 - 3.2983e+17j) * hz_scale)
]
Au_IIR2 = Medium(epsilon=1, E_susceptibilities=Au_IIR2_susc)

#------------------------------------------------------------------
# gold (Au) square root IIR 500-1000 nm model
Au_sqrt_IIR2_susc = [
	IIR_Susceptibility((-8.6705e+15) * hz_scale, (-4.3885e+15) * hz_scale),
	IIR_Susceptibility((-8.0645e+14 + 4.4853e+15j) * hz_scale, (8.3663e+14 + 2.3093e+14j) * hz_scale),
	IIR_Susceptibility((-6.3678e+13) * hz_scale, (6.1786e+15) * hz_scale)
]
Au_sqrt_IIR2 = Medium(epsilon=1, E_susceptibilities=Au_sqrt_IIR2_susc)

#------------------------------------------------------------------
# fictious 1
fic1_frq0 = 3e14 * hz_scale
fic1_gam0 = 0.5e14 * hz_scale
fic1_sig0 = 3

fic1_frq1 = 5e14 * hz_scale
fic1_gam1 = 1e14 * hz_scale
fic1_sig1 = 3
fic_susc = [LorentzianSusceptibility(frequency=fic1_frq0, gamma=fic1_gam0, sigma=fic1_sig0),
            LorentzianSusceptibility(frequency=fic1_frq1, gamma=fic1_gam1, sigma=fic1_sig1)]

fic = Medium(epsilon=1, mu=1, E_susceptibilities=fic_susc)

#------------------------------------------------------------------
# elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83, 1998
# wavelength range: 0.2 - 12.4 um

# metal_range = mp.FreqRange(min=um_scale/12.4, max=um_scale/0.2)

# silver (Ag)

Ag_plasma_frq = 9.01*eV_um_scale
Ag_f0 = 0.845
Ag_frq0 = 1e-10
Ag_gam0 = 0.048*eV_um_scale
Ag_sig0 = Ag_f0*Ag_plasma_frq**2/Ag_frq0**2
Ag_f1 = 0.065
Ag_frq1 = 0.816*eV_um_scale      # 1.519 um
Ag_gam1 = 3.886*eV_um_scale
Ag_sig1 = Ag_f1*Ag_plasma_frq**2/Ag_frq1**2
Ag_f2 = 0.124
Ag_frq2 = 4.481*eV_um_scale      # 0.273 um
Ag_gam2 = 0.452*eV_um_scale
Ag_sig2 = Ag_f2*Ag_plasma_frq**2/Ag_frq2**2
Ag_f3 = 0.011
Ag_frq3 = 8.185*eV_um_scale      # 0.152 um
Ag_gam3 = 0.065*eV_um_scale
Ag_sig3 = Ag_f3*Ag_plasma_frq**2/Ag_frq3**2
Ag_f4 = 0.840
Ag_frq4 = 9.083*eV_um_scale      # 0.137 um
Ag_gam4 = 0.916*eV_um_scale
Ag_sig4 = Ag_f4*Ag_plasma_frq**2/Ag_frq4**2

Ag_susc = [DrudeSusceptibility(frequency=Ag_frq0, gamma=Ag_gam0, sigma=Ag_sig0),
           LorentzianSusceptibility(frequency=Ag_frq1, gamma=Ag_gam1, sigma=Ag_sig1),
           LorentzianSusceptibility(frequency=Ag_frq2, gamma=Ag_gam2, sigma=Ag_sig2),
           LorentzianSusceptibility(frequency=Ag_frq3, gamma=Ag_gam3, sigma=Ag_sig3),
           LorentzianSusceptibility(frequency=Ag_frq4, gamma=Ag_gam4, sigma=Ag_sig4)]

Ag = Medium(epsilon=1.0, E_susceptibilities=Ag_susc)

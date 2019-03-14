import ffip

# default unit length is 1 um
um_scale = 1.0

# conversion factor for eV to 1/um [=1/hc]
eV_um_scale = um_scale/1.23984193

#------------------------------------------------------------------
# gold (Au)

Au_plasma_frq = 9.03*eV_um_scale
Au_f0 = 0.760
Au_frq0 = 1e-10
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

Au_susc = [ffip.DrudeSusceptibility(frequency=Au_frq0, gamma=Au_gam0, sigma=Au_sig0),
           ffip.LorentzianSusceptibility(frequency=Au_frq1, gamma=Au_gam1, sigma=Au_sig1),
           ffip.LorentzianSusceptibility(frequency=Au_frq2, gamma=Au_gam2, sigma=Au_sig2),
           ffip.LorentzianSusceptibility(frequency=Au_frq3, gamma=Au_gam3, sigma=Au_sig3),
           ffip.LorentzianSusceptibility(frequency=Au_frq4, gamma=Au_gam4, sigma=Au_sig4)]

Au = ffip.Medium(epsilon=1.0, E_susceptibilities=Au_susc)
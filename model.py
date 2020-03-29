#
# Let us try to reproduce the fit from Allen et al (2003)
#
# The LF is given, the model is extremely simple, M-L relation without any scatter.
#
# The best fit cosmology is:
#
# sigma8 = 0.510 \pm 0.019 Omegam^(-0.253\pm0.024)
#
# Since the measured LF is at LCDM with Omegam=0.3, we will have to use a value
# of sigma8 commensurate with the above best fitin order to get a decent fit,
# and make sure that we do not have change the LF or the mass-luminosity
# relation.
# 
# Thus for Omegam=0.3, sigma8 = 0.6916
#
# The mass-luminosity relation is:
# 
# log10(E(z)M200/(h50inv Msun)) = alpha * log10(L/E(z)/10^{44}/h50inv^2/erg/s^-1) + log10(M0/(h50inv Msun))

#
# This is not the best fit cosmology without the scatter incorporated, but let us forget about it for the moment.
#
from math import log10
import sys
sys.path.append("../aum/aum/divyacheck/lib64/python3.7/site-packages/")
import cosmology as cc
import pandas
import numpy as np

# Allen et al: Use 0.6916 as sigma8, Use T0=2.726, h=0.72, ns=1.0, Omegab = 0.03858
c_allen = cc.cosmology(0.3,0.0,-1.0,0.0,0.03858,0.72,2.726,0.6916,1.0,log10(8.0),1.0)

z_allen = 0.21

#
# dn/dL = dn/dM * dM/dL
# M200 = M0 L^{alpha}/E(z)^{alpha+1}
# dM/dL = \alpha M200/L
#
alpha = 0.76
logM0 = 14.29

#
# Now compute the luminosity function with the above simple formula
#
Lxarr = np.logspace(np.log10(10.), np.log10(30.0), 100)
Ez = c_allen.Eofz(z_allen)
Marr = 10.**(logM0 + alpha*(np.log10(Lxarr/Ez)) - np.log10(Ez))

# The above mass is in units of h50inv Msun and Lx in units of 10^44 h50^{-2} erg s^{-1}
hval = 0.50

# Get masses in hinv Msun so that they can be passed to my code
Marrhinv = Marr*0.50

Phiarr = Lxarr * 0.0
for ii in range(Lxarr.size):
    # Get the mass function from the code in hinv^{-4} Mpc^{-3} Msun^{-1} units
    # and convert to appropriate h50 units mass function
    mf_h50 = c_allen.MF_Jenkins(Marrhinv[ii], z_allen) * hval**4

    # Now multiply by dM/dL all in h50 related units
    Phiarr[ii] = mf_h50*alpha*Marr[ii]/Lxarr[ii]

lum, nclus, lf, elf = np.loadtxt("Allen_etal_lf.dat", unpack=True)

import pylab as pl

ax = pl.subplot(111)
ax.errorbar(lum, lf, elf, label="Allen et al. (2003)", fmt="b.")
ax.plot(Lxarr, Phiarr, label="Model")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$L_X$ ($10^{44} h_{50}^{-2}{\rm erg}{\rm s^{-1}}$)")
ax.set_ylabel(r"$\Phi(L_X)$ $(h_{50}^{5} {\rm Mpc}^{-3}/{\rm 10^{44} erg s^{-1}})$")
pl.savefig("LF.pdf")

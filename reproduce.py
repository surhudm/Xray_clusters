#
# Let us try to reproduce the numbers of BCS clusters quoted in Allen et al. (2003)
#
# The data is from Ebeling et al. (1998, 2000)
#
# The cosmology assumed for the computed Lx in the data file is Omegam=1, q0=0.5, OmegaL=0.0 (see references)
#
from math import log10
import sys
sys.path.append("../aum/aum/divyacheck/lib64/python3.7/site-packages/")
import cosmology as cc
import pandas
import numpy as np

# Ebeling
c_ebeling = cc.cosmology(1.0,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)
c_allen = cc.cosmology(0.3,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0)

# Read in data file
df_ebeling = pandas.read_csv("clusters.dat", delim_whitespace=1)
df_ebeling_e = pandas.read_csv("e_clusters.dat", delim_whitespace=1)
df_ebeling["Lx_allen"] = df_ebeling["lx(44)"].values*0.0
df_ebeling_e["Lx_allen"] = df_ebeling_e["lx(44)"].values*0.0

for ii in range(df_ebeling.Lx_allen.values.size):
    # Find ratio between Dl in Allen versus Ebeling et al. cosmology irrespective of H0 values (that is what aum gives us)
    # Ebeling: 50 km/s/Mpc
    # Allen: 50 km/s/Mpc
    ratio = c_allen.Dlofz(df_ebeling.z.values[ii])**2/c_ebeling.Dlofz(df_ebeling.z.values[ii])**2

    df_ebeling["Lx_allen"].values[ii] = df_ebeling["lx(44)"].values[ii]*ratio

for ii in range(df_ebeling_e.Lx_allen.values.size):
    # Find ratio between Dl in Allen versus Ebeling et al. cosmology irrespective of H0 values (that is what aum gives us)
    # Ebeling: 50 km/s/Mpc
    # Allen: 50 km/s/Mpc
    ratio = c_allen.Dlofz(df_ebeling_e.z.values[ii])**2/c_ebeling.Dlofz(df_ebeling_e.z.values[ii])**2

    df_ebeling_e["Lx_allen"].values[ii] = df_ebeling_e["lx(44)"].values[ii]*ratio

# Ok now make the cuts and figure out the numbers
def return_sample_idx(Lav, Lplus, Lminus):
    idx_e = (df_ebeling_e.z.values<=0.3) & (df_ebeling_e["Lx_allen"].values>Lav-Lminus) & (df_ebeling_e["Lx_allen"].values<=Lav+Lplus)
    idx = (df_ebeling.z.values<=0.3) & (df_ebeling["Lx_allen"].values>Lav-Lminus) & (df_ebeling["Lx_allen"].values<=Lav+Lplus)
    Nclusters = np.sum(idx) + np.sum(idx_e)
    print ("Number of clusters:", Nclusters)
    Lxclusters = np.zeros(Nclusters)
    Lxclusters[:np.sum(idx)] = df_ebeling["Lx_allen"].values[idx]
    Lxclusters[np.sum(idx):] = df_ebeling_e["Lx_allen"].values[idx_e]
    print ("Mean Lx:", np.mean(Lxclusters))
    zclusters = np.zeros(Nclusters)
    zclusters[:np.sum(idx)] = df_ebeling["z"].values[idx]
    zclusters[np.sum(idx):] = df_ebeling_e["z"].values[idx_e]
    print ("Mean redshift:", np.mean(zclusters))
    return idx, idx_e, Lxclusters, zclusters

idx, idx_e, Lxclusters, zclusters = return_sample_idx(11.73, 99.0, 1.73)

# Plot the clusters to see where the problematic cluster may lie
import pylab as pl

ax = pl.subplot(111)
ax.scatter(df_ebeling["z"].values, df_ebeling["Lx_allen"].values)
ax.scatter(df_ebeling_e["z"].values, df_ebeling_e["Lx_allen"].values)
ax.scatter(df_ebeling["z"].values[idx], df_ebeling["Lx_allen"].values[idx])
ax.scatter(df_ebeling_e["z"].values[idx_e], df_ebeling_e["Lx_allen"].values[idx_e])
ax.set_yscale("log")
ax.axvline(0.3)
ax.set_ylim(9.0)

# From the plot it is clear that at the first two BCS luminosity bin, the flux
# limit clearly plays a role. Figure out the flux limit and work backwards to
# obtain the Luminosity-max redshift curve.

Flim = 2.8e-12
zarr = np.linspace(0.01, 1.0, 100)
Larr = zarr*0.0
Varr = zarr*0.0
Mpc_to_cm = 3.08567758128E+24
for ii in range(zarr.size):
    Larr[ii] = 4.0*np.pi*Flim*Mpc_to_cm**2*c_allen.Dlofz(zarr[ii])**2.0/0.50**2/1e44
    Varr[ii] = 4.14/3*c_allen.Dcofz(zarr[ii])**3/0.50**3

maxVmax = 4.14/3*c_allen.Dcofz(0.3)**3/0.50**3

ax.plot(zarr, Larr)

pl.savefig("Allen_2003.png")
pl.clf()

# This plot looks reasonable, there are few oddities, but my guess is they are
# all related to details about k-corrections, etc. In what follows we ignore
# those details for the time being.

lav, lplus, lminus, Num, Phi, Phi_err = np.loadtxt("Allen_etal_lf.dat", unpack=True)
bcslav = lav[0:3]
bcslplus = lplus[0:3]
bcslminus = lminus[0:3]
bcsNum = Num[0:3]
bcsPhi = Phi[0:3]
bcsPhi_err = Phi_err[0:3]

# Spline for zmax to Larr
from scipy.interpolate import interp1d
spline = interp1d(Larr, Varr)

# Now for every cluster figure out the Vmax
Vmax = spline(Lxclusters)
cutidx = Vmax>maxVmax
Vmax[cutidx]=maxVmax

# Get the luminosity function and compare with Allen et al.
bins = bcslav-bcslminus
bins = np.append(bins, bcslav[2]+bcslplus[2])
print(bins)

hist, bin_edges = np.histogram(Lxclusters, bins, weights=1/Vmax)
xhist = hist*0.0
for ii in range(hist.size):
    cutidx = (Lxclusters>=bcslav[ii]-bcslminus[ii]) & (Lxclusters<bcslav[ii]+bcslplus[ii])
    print (np.sum(cutidx))
    xhist[ii] = np.average(Lxclusters[cutidx])

# Now divide by Delta L to get Phi(L)
hist = hist/(bcslplus+bcslminus)

ax = pl.subplot(111)
ax.errorbar(bcslav, bcsPhi, bcsPhi_err, fmt=".", label="BCS Allen et al.")
ax.scatter(xhist, hist, label="BCS my estimate")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlim(8.0, 30.0)
ax.set_ylim(1e-11, 1e-8)
ax.set_xlabel(r"$L_{\rm x}$ ($10^{44}h_{50}^{-2}$ erg$s^{-1}$)")
ax.set_ylabel("$\Phi$ ($h_{50}^{-5}$ Mpc$^{-3}$ /$10^{44}$s/erg)")
ax.legend()
pl.savefig("Phi.png")

df_ebeling.to_csv("clusters_wallen.dat", sep=" ", index=False)
df_ebeling_e.to_csv("e_clusters_wallen.dat", sep=" ", index=False)

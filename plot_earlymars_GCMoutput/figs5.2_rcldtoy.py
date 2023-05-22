# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font',family = 'serif',size = 12)

from scipy.io import netcdf
from scipy.interpolate import interp1d

rp  = 3389.5e3 #6371e3
rcp = 8.31/32e-3 / 1062.5
spin= 7.292e-5 
grav= 3.71 #9.8
ro2 = 8.31/32e-3 
Lv = 2.5e6
Rv = 461.5

def esat(T):
    return 610.78 * np.exp(-Lv/Rv*(1./T - 1./273.16))
    
def area_mean(var2d, lat): 
    clat = np.cos(lat / 180.*np.pi)
    var1 = np.mean(var2d, axis=1) *clat 
    var2 = np.sum(var1, axis=0) / np.nansum(clat)
    return var2

def vertical_interp(var3d, p3d, p1d): #pres, lat, lon, for visualization
    aa, nlat, nlon = var3d.shape
    varout = np.zeros((len(p1d),nlat,nlon))
    for i in range(nlat):
        for j in range(nlon): 
            ff = interp1d(np.log(p3d[:,i,j]), var3d[:,i,j], 
                    bounds_error=False, fill_value=np.nan)
            varout[:,i,j] = ff(np.log(p1d))
    return varout

def vertical_integration(var3d, phalf3d):
    npfull, nlat, nlon = var3d.shape
    varout = np.zeros((nlat, nlon))
    for i in range(npfull):
        varout[:,:] = varout[:,:] + var3d[i,:,:]* \
            (phalf3d[i+1,:,:] - phalf3d[i,:,:])/grav
    return varout

pk = np.array([219.4067, 489.5209, 988.2418, 1805.201, 2983.724, 4462.334, 6160.587, 
    7851.243, 7731.271, 7590.131, 7424.086, 7228.744, 6998.933, 6728.574, 
    6410.509, 6036.322, 5596.111, 5078.225, 4468.96, 3752.191, 2908.949, 
    2084.739, 1334.443, 708.499, 252.136, 0, 0])
bk = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0.01505309, 0.03276228, 0.05359622, 0.07810627, 
    0.1069411, 0.1408637, 0.180772, 0.227722, 0.2829562, 0.3479364, 
    0.4243822, 0.5143168, 0.6201202, 0.7235355, 0.8176768, 0.8962153, 
    0.9534761, 0.9851122, 1])

def evap_flux(ts):
    ta = ts/2**0.25
    farea = 7.596e-3
    return 1e-3 *(esat(ts)) / (Rv * ts) *farea

def iplus(ts, taui, pc): #pressure is nondimensional
#     pc = 0.2
    tc = ts*pc**rcp
    pres = np.logspace(np.log10(pc),0,1000)
    lnp  = np.log(pres)
    iplus_list = []

#     for taui in tauinf:
#         trans= np.exp(-taui*(pres-pc))
    trans= np.exp(-taui*(pres**2-pc**2))
    tair = ts*pres**rcp
    integrand = trans*4*5.67e-8*tair**4 *rcp
    iplus= 5.67e-8*tc**4 + np.sum((integrand[:-1]+integrand[1:])/2.*(lnp[1:]-lnp[:-1]))
    return iplus


fs = 109.9644*(1-0.31) 
tauco2 = 4.
pc = 0.2
# e2 = 0.9 #co2
# taucon_list = np.array([9e2, 9e3, 9e4, 9e5, 2e6, 3e6, 9e6])
taucon_list = np.logspace(5,7,100)
# ts_list     = []

print (fs/5.67e-8)**0.25

# mass_to_depth = 3./4/(1e3*10e-6)
rcld_list = np.array([1, 5, 10, 20])*1e-6
# absmean_list = np.array([0.12365507537923022, 0.57890201796958596, 
#     0.85502030303407794, 1.0509896784057426]) #Planck mean
absmean_list = np.array([0.056257576102087141, 0.36965137293495259, 
    0.71150878241861992, 0.9908548548054491]) #Rosseland mean
mass_to_depth = absmean_list *3 /4 /(1e3*rcld_list)
print mass_to_depth

labels = [r'r=1$\mu$m',r'r=5$\mu$m',r'r=10$\mu$m',r'r=20$\mu$m']
colors = ['g', 'b','r','k']
fig = plt.figure(figsize=(4,3))
# cmap = mpl.cm.RdBu_r 
ax = fig.add_subplot(111)

for i in range(len(mass_to_depth)):
    ts_list = []
    for taucon in taucon_list:
        # taucon = 1e3
        ts0 = 220.
        dts = 1.
        while abs(dts)>1e-2:

            tc = ts0*pc**rcp
            mcld = evap_flux(ts0) * taucon
            e1 = 1.- np.exp(-mass_to_depth[i]*mcld)
        #         tsurf = (fs/5.67e-8 *(4.-e1*e2)/(2-e1)/(2-e2))**0.25
        #         tsurf = (fs*2./(2-emis)/5.67e-8)**0.25
            fco2 = iplus(ts0, tauco2, pc)
            olr  = fco2*(1-e1)+e1*5.67e-8 *tc**4
    #         print olr, fs, ts0

            dts = fs - olr
            ts0 = ts0 + dts
        ts_list.append(ts0)

    line1, = ax.semilogx(taucon_list/86400, ts_list, color=colors[i],
        label=labels[i])

ax.set_xlabel(r'Cloud lifetime $\tau_c$ (day)')
ax.set_ylabel(r'Surface temperature (K)')
# ax.set_ylabel(r'Global mean cloud path (kg m$^{-2}$)')
ax.grid()
ax.set_xlim([1,100])
ax.legend(loc=0,fontsize=10)
# ax.legend(numpoints=1,frameon=False,
#     bbox_to_anchor=(0.1, -0.4, 1.3, .102), loc='lower left',
#     ncol=2, mode="expand", borderaxespad=0.)
ax.tick_params(which='both',direction='out')

# ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
# color = 'blue'
# ax2.set_ylabel('Cloud optical thickness (COT)', color=color)  # we already handled the x-label with ax1
# line2, = ax2.loglog(taucon_list/86400, evap_flux(np.array(ts_list)) * taucon_list*mass_to_depth, \
#     '--',color=color, label='cloud optical depth')
# ax2.plot([1e-2, 2e2], [1,1], color='y', linewidth=5, alpha=0.5)
# ax2.tick_params(axis='y', labelcolor=color)
# # ax.set_ylim([1e-5, ])
# # ax2.legend(numpoints=1,frameon=False,
# #     bbox_to_anchor=(0.1, -0.5, 1.3, .102), loc='lower left',
# #     ncol=2, mode="expand", borderaxespad=0.)
# ax2.legend([line1, line2], ['Surface temperature', 'Cloud optical thickness'], 
#     numpoints=1, loc='upper left', fontsize=10)
# ax2.set_xlim([1e-2,2e2])
# ax2.set_ylim([1e-5,3e2])
# ax2.tick_params(which='both',direction='out')
# ax.text(0, 1.05, '(b) Toy model',transform=ax.transAxes, fontsize=12)
# ax.text(0.82, 0.58, 'COT=1', transform=ax.transAxes,
#     color='y',fontsize=10)

plt.show()
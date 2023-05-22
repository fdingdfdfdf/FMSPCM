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

tau_auto = np.array([9e2, 9e3, 9e4, 9e5, 1.5e6, 3e6, 9e6])
sflux_mean  = tau_auto *0.
rain_mean   = tau_auto *0.
cld_mean    = tau_auto *0.
dirname  = ['1_sh9e2/day4122h00.',
    '2_sh9e3/day5496h00.',
    '3_sh9e4/day2748h00.',
    '4_sh9e5/day5496h00.',
    '5_sh1.5e6/day8244h00.',
    '6_sh3e6/day17175h00.',
    '7_sh9e6/day12366h00.']
for i in range(len(dirname)):
    filename = '../data/SH/'+dirname[i]+'atmos_avg.nc'
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    sflux  = f.variables['flux_lhe'].data[0,:,:]
    rain   = f.variables['total_rain'].data[0,:,:]
    cloud = f.variables['cld_amt'].data[0,:,:,:] #pres, lat, lon

    sflux_mean[i] = area_mean(sflux, lat)
    rain_mean[i]  = area_mean(rain, lat)

    ps    = f.variables['ps'].data[0,:,:]
    phalf = np.zeros((len(pk),len(lat),len(lon)))
    pres  = cloud *1.
    for ij in range(len(pk)):
        phalf[ij,:,:] = ps[:,:]*bk[ij] + pk[ij]
    pres[:,:,:] = (phalf[1:,:,:] - phalf[:-1,:,:])  \
        / (np.log(phalf[1:,:,:]) - np.log(phalf[:-1,:,:]))
    cld_mean[i] = area_mean(vertical_integration(cloud, phalf), lat)

    f.close()
mass_to_depth = 3./4*1e2 #from cloud mass to optical thickness
fig = plt.figure(figsize=(8.5,3))
cmap = mpl.cm.RdBu_r 
ax = fig.add_subplot(121)
ax.loglog(tau_auto/86400, sflux_mean/Lv*tau_auto* mass_to_depth,
    'k',linestyle='none', marker='s',label='Evaporation, SP')
# ax.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
#     '--',label='Precipitation')
ax.loglog(tau_auto/86400, cld_mean* mass_to_depth, 'k-', 
    label='Cloud mass, SP')


print sflux_mean, rain_mean

tau_auto = np.array([9e2, 9e3, 9e4, 9e5,1.5e6, 3e6, 9e6])
sflux_mean  = tau_auto *0.
rain_mean   = tau_auto *0.
cld_mean    = tau_auto *0.
dirname  = ['1_eq9e2/day5496h00.',
            '2_eq9e3/day5496h00.',
            '3_eq9e4/day2748h00.',
            '4_eq9e5/day5496h00.',
            '5_eq1.5e6/day10992h00.',
            '6_eq3e6/day10305h00.',
            '7_eq9e6/day7557h00.']
ts_mean  = tau_auto *0.
for i in range(len(dirname)):
    filename = '../data/EQ/'+dirname[i]+'atmos_avg.nc'       
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    sflux  = f.variables['flux_lhe'].data[0,:,:]
    rain   = f.variables['total_rain'].data[0,:,:]
    cloud = f.variables['cld_amt'].data[0,:,:,:] #pres, lat, lon

    sflux_mean[i] = area_mean(sflux, lat)
    rain_mean[i]  = area_mean(rain, lat)

    ps    = f.variables['ps'].data[0,:,:]
    phalf = np.zeros((len(pk),len(lat),len(lon)))
    pres  = cloud *1.
    for ij in range(len(pk)):
        phalf[ij,:,:] = ps[:,:]*bk[ij] + pk[ij]
    pres[:,:,:] = (phalf[1:,:,:] - phalf[:-1,:,:])  \
        / (np.log(phalf[1:,:,:]) - np.log(phalf[:-1,:,:]))
    cld_mean[i] = area_mean(vertical_integration(cloud, phalf), lat)

    f.close()
mass_to_depth = 3./4*1e2 #from cloud mass to optical thickness

ax.loglog(tau_auto/86400, sflux_mean/Lv*tau_auto* mass_to_depth,
    'r',linestyle='none',marker='o',label='Evaporation, EQ')
# ax.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
#     '--',label='Precipitation')
ax.loglog(tau_auto/86400, cld_mean* mass_to_depth, 'r--', 
    label='Cloud mass, EQ')
ax.plot([1e-2, 2e2], [1,1], color='y', linewidth=5, alpha=0.5)
ax.plot([45./24/60, 45./24/60], [1e-5,1e2], color='b', linewidth=5, alpha=0.5)
ax.set_xlabel(r'Cloud lifetime $\tau_c$ (day)')
#ax.set_ylabel(r'Global mean cloud path (kg m$^{-2}$)')
ax.set_ylabel('Global mean \n cloud optical thickness (COT)')
ax.grid()
# ax.legend(numpoints=1,frameon=False,
#     bbox_to_anchor=(-0.3, -0.5, 1.5, .102), loc='lower left',
#     ncol=2, mode="expand", borderaxespad=0.)
ax.legend(numpoints=1, loc='upper left', fontsize=10)
ax.set_xlim([1e-2,2e2])
ax.tick_params(which='both',direction='out')
ax.text(0, 1.05, '(a) GCM',transform=ax.transAxes, fontsize=12)
ax.text(0.82, 0.63, 'COT=1', transform=ax.transAxes,
    color='y',fontsize=10)
ax.text(0.15, 0.1, 'Earth', transform=ax.transAxes,
    color='b',fontsize=10)
    
print sflux_mean, rain_mean

# plt.tight_layout()
# fig.savefig('SFLUX_CLD_EQ.pdf')


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
taucon_list = np.logspace(3,7,100)
ts_list     = []

print (fs/5.67e-8)**0.25
for taucon in taucon_list:
    # taucon = 1e3
    ts0 = 220.
    dts = 1.
    while abs(dts)>1e-2:

        tc = ts0*pc**rcp
        mcld = evap_flux(ts0) * taucon
        e1 = 1.- np.exp(-75.*mcld)
    #         tsurf = (fs/5.67e-8 *(4.-e1*e2)/(2-e1)/(2-e2))**0.25
    #         tsurf = (fs*2./(2-emis)/5.67e-8)**0.25
        fco2 = iplus(ts0, tauco2, pc)
        olr  = fco2*(1-e1)+e1*5.67e-8 *tc**4
#         print olr, fs, ts0

        dts = fs - olr
        ts0 = ts0 + dts
    ts_list.append(ts0)

# fig = plt.figure(figsize=(4,3))
# cmap = mpl.cm.RdBu_r 
ax = fig.add_subplot(122)
line1, = ax.semilogx(taucon_list/86400, ts_list,
    'k-',label='Surface temperature')

ax.set_xlabel(r'Cloud lifetime $\tau_c$ (day)')
ax.set_ylabel(r'Surface temperature (K)')
# ax.set_ylabel(r'Global mean cloud path (kg m$^{-2}$)')
ax.grid()
# ax.legend(numpoints=1,frameon=False,
#     bbox_to_anchor=(0.1, -0.4, 1.3, .102), loc='lower left',
#     ncol=2, mode="expand", borderaxespad=0.)
ax.tick_params(which='both',direction='out')

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
color = 'blue'
ax2.set_ylabel('Cloud optical thickness (COT)', color=color)  # we already handled the x-label with ax1
line2, = ax2.loglog(taucon_list/86400, evap_flux(np.array(ts_list)) * taucon_list*mass_to_depth, \
    '--',color=color, label='cloud optical depth')
ax2.plot([1e-2, 2e2], [1,1], color='y', linewidth=5, alpha=0.5)
ax2.tick_params(axis='y', labelcolor=color)
# ax.set_ylim([1e-5, ])
# ax2.legend(numpoints=1,frameon=False,
#     bbox_to_anchor=(0.1, -0.5, 1.3, .102), loc='lower left',
#     ncol=2, mode="expand", borderaxespad=0.)
ax2.legend([line1, line2], ['Surface temperature', 'Cloud optical thickness'], 
    numpoints=1, loc='upper left', fontsize=10)
ax2.set_xlim([1e-2,2e2])
ax2.set_ylim([1e-5,3e2])
ax2.tick_params(which='both',direction='out')
ax.text(0, 1.05, '(b) Toy model',transform=ax.transAxes, fontsize=12)
ax.text(0.82, 0.58, 'COT=1', transform=ax.transAxes,
    color='y',fontsize=10)

plt.show()
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

fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True)
# fig.set_size_inches(12, 6)
cmap = mpl.cm.spectral #RdBu_r 
levelp = np.logspace(2, np.log10(1e5), 40)

titles = [r'(a) SP, $\tau_c$=0.01 day',r'(b) SP, $\tau_c$=10 day',r'(c) SP, $\tau_c$=100 day']
tau_auto = np.array([9e2, 9e5, 9e6])
sflux_mean  = tau_auto *0.
rain_mean   = tau_auto *0.
cld_mean    = tau_auto *0.
dirname  = ['1_sh9e2/day4122h00.',
    '4_sh9e5/day5496h00.',
    '7_sh9e6/day12366h00.']
for i in range(len(dirname)):
    filename = '/n/wordsworth_lab/fding/fms_output/fv_mars_shcloudevap'+dirname[i]+'atmos_avg.nc'
    filename = '../data/SH/'+dirname[i]+'atmos_avg.nc'
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
#     sflux  = f.variables['flux_lhe'].data[0,:,:]
#     rain   = f.variables['total_rain'].data[0,:,:]
    cloud = f.variables['cld_amt'].data[0,:,:,:] #pres, lat, lon
#     sflux_mean[i] = area_mean(sflux, lat)
#     rain_mean[i]  = area_mean(rain, lat)
    ps    = f.variables['ps'].data[0,:,:]
    phalf = np.zeros((len(pk),len(lat),len(lon)))
    pres  = cloud *1.
    for ij in range(len(pk)):
        phalf[ij,:,:] = ps[:,:]*bk[ij] + pk[ij]
    pres[:,:,:] = (phalf[1:,:,:] - phalf[:-1,:,:])  \
        / (np.log(phalf[1:,:,:]) - np.log(phalf[:-1,:,:]))
#     cld_mean[i] = area_mean(vertical_integration(cloud, phalf), lat)

#plot
    X,Y = np.meshgrid(lat, levelp)
    y1out  = vertical_interp(cloud, pres, levelp)
    y1out[y1out<2e-15] = 2e-15
    ax = axes.flat[i]
    CS2 = ax.contourf(X,Y, np.mean(y1out,axis=2), locator=mpl.ticker.LogLocator(), 
        levels=np.logspace(-15,-4,22+1), cmap=cmap  )
    ax.set_yscale('log')
    ax.set_ylim([1e5, 4e2])
    if (i==0):
        ax.set_ylabel(r'Pressure (Pa)')
    ax.tick_params(which='both',direction='out')
    ax.grid()
    
    ax.text(-90, 3e2, titles[i], fontsize=12)
    f.close()
# mass_to_depth = 3./4*1e2 #from cloud mass to optical thickness

titles = [r'(d) EQ, $\tau_c$=0.01 day',r'(e) EQ, $\tau_c$=10 day',r'(f) EQ, $\tau_c$=100 day']
dirname  = ['1_eq9e2/day5496h00.',
            '4_eq9e5/day5496h00.',
            '7_eq9e6/day7557h00.']
sflux_mean  = tau_auto *0.
rain_mean   = tau_auto *0.
cld_mean    = tau_auto *0.
for i in range(len(dirname)):
    filename = '../data/EQ/'+dirname[i]+'atmos_avg.nc'   
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
#     sflux  = f.variables['flux_lhe'].data[0,:,:]
#     rain   = f.variables['total_rain'].data[0,:,:]
    cloud = f.variables['cld_amt'].data[0,:,:,:] #pres, lat, lon
#     sflux_mean[i] = area_mean(sflux, lat)
#     rain_mean[i]  = area_mean(rain, lat)
    ps    = f.variables['ps'].data[0,:,:]
    phalf = np.zeros((len(pk),len(lat),len(lon)))
    pres  = cloud *1.
    for ij in range(len(pk)):
        phalf[ij,:,:] = ps[:,:]*bk[ij] + pk[ij]
    pres[:,:,:] = (phalf[1:,:,:] - phalf[:-1,:,:])  \
        / (np.log(phalf[1:,:,:]) - np.log(phalf[:-1,:,:]))
#     cld_mean[i] = area_mean(vertical_integration(cloud, phalf), lat)

#plot
    X,Y = np.meshgrid(lat, levelp)
    y1out  = vertical_interp(cloud, pres, levelp)
    y1out[y1out<2e-15] = 2e-15
    ax = axes.flat[i+3]
    CS2 = ax.contourf(X,Y, np.mean(y1out,axis=2), locator=mpl.ticker.LogLocator(), 
        levels=np.logspace(-15,-4,22+1), cmap=cmap  )
    ax.set_yscale('log')
    ax.set_ylim([1e5, 4e2])
    ax.set_xlabel('Lat (deg)')
    if (i==0):
        ax.set_ylabel(r'Pressure (Pa)')
    ax.grid()
    ax.tick_params(which='both',direction='out')
    ax.text(-90, 3e2, titles[i], fontsize=12)
    f.close()

CB = fig.colorbar(CS2, ax=axes.ravel().tolist(), pad=0.1, location='bottom')
CB.ax.set_xlabel(r'Water cloud mass concentration (kg/kg)')
    # CB.set_ticks()
    
# plt.tight_layout()
# fig.savefig('CLD_MASS_EQ.pdf')

plt.show()
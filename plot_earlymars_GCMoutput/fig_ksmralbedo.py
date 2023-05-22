# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font',family = 'serif',size = 12)

from scipy.io import netcdf
from scipy.interpolate import interp1d

rp  = 3389.5e3 #6371e3
# rcp = 8.31/32e-3 / 1062.5
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

# pk = np.array([219.4067, 489.5209, 988.2418, 1805.201, 2983.724, 4462.334, 6160.587, 
#     7851.243, 7731.271, 7590.131, 7424.086, 7228.744, 6998.933, 6728.574, 
#     6410.509, 6036.322, 5596.111, 5078.225, 4468.96, 3752.191, 2908.949, 
#     2084.739, 1334.443, 708.499, 252.136, 0, 0])
# bk = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0.01505309, 0.03276228, 0.05359622, 0.07810627, 
#     0.1069411, 0.1408637, 0.180772, 0.227722, 0.2829562, 0.3479364, 
#     0.4243822, 0.5143168, 0.6201202, 0.7235355, 0.8176768, 0.8962153, 
#     0.9534761, 0.9851122, 1])
mass_to_depth = 3./4 /(1e3*10e-6) #from cloud mass to optical thickness

fig = plt.figure(figsize=(4,5))
cmap = mpl.cm.RdBu_r 

ax = fig.add_subplot(311)
ax2= fig.add_subplot(212) 

#KSMR21, more ice
filename = '../data/KITE21/8_kite4e6iceplus/day2748h00.atmos_static.nc'   
f = netcdf.netcdf_file(filename, 'r')
lat = f.variables['grid_yt'].data
lon = f.variables['grid_xt'].data
zsurf  = f.variables['zsurf'].data[:,:]
f.close()

filename = '../data/KITE21/8_kite4e6iceplus/day2748h00.atmos_echeck.nc'   
f = netcdf.netcdf_file(filename, 'r')
# lat = f.variables['grid_yt'].data
# lon = f.variables['grid_xt'].data
solar  = f.variables['solar'].data[0,:,:]
osr  = f.variables['OSR'].data[0,:,:]
f.close()

water = zsurf*0.
for i in range(len(lat)):
    for j in range(len(lon)):
        if (zsurf[i,j]>1e3 and (lat[i])<60):
            water[i,j] = 1.

X,Y = np.meshgrid(lon, lat)
CS2 = ax.contour(X,Y,zsurf/1e3,colors='k',linewidths=0.4)
ax.contourf(X,Y,water,levels=[0,0.6,1],colors=['sandybrown','w'])#,zorder=2.5)
ax.tick_params(which='both',direction='out')
ax.set_ylabel('Lat (deg)')
# ax.set_ylabel('Latitude (deg)')
# ax.set_xlabel('Longitude (deg)')
# ax.text(0, 1.05, '(b) EQ', fontsize=12, transform=ax.transAxes)
CS2 = ax2.contourf(X,Y,osr/solar,cmap = cmap)
ax2.tick_params(which='both',direction='out')
ax2.set_ylabel('Lat (deg)')
ax2.set_xlabel('Lon (deg)')

CB = fig.colorbar(CS2, ax=ax2, pad=0.3, orientation='horizontal',
    ticks = np.linspace(0.8,0.95,4))
CB.ax.set_xlabel(r'Planetary albedo')
# # plt.tight_layout()
# # fig.savefig('TS_TAU_EQ.pdf')


fig = plt.figure(figsize=(4,3))
cmap = mpl.cm.RdBu

ax = fig.add_subplot(111)

#KSMR21, more ice
filename = '../data/KITE21/8_kite4e6iceplus/day2748h00.atmos_avg.nc'   
f = netcdf.netcdf_file(filename, 'r')
lat = f.variables['grid_yt'].data
lon = f.variables['grid_xt'].data
evap  = f.variables['flux_lhe'].data[0,:,:]
rain  = f.variables['total_rain'].data[0,:,:]
f.close()

X,Y = np.meshgrid(lon, lat)
CS2 = ax.contourf(X,Y,rain-evap,cmap = cmap,
    norm=mpl.colors.Normalize(vmin=-4, vmax=4))
ax.tick_params(which='both',direction='out')
ax.set_ylabel('Lat (deg)')
ax.set_xlabel('Lon (deg)')

CB = fig.colorbar(CS2, ax=ax, pad=0.3, orientation='horizontal')
CB.ax.set_xlabel(r'Precipitation-Evaporation (W m$^{-2}$)')

plt.show()
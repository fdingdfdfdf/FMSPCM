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

# tau_auto = np.array([9e2, 9e3, 9e4, 9e5, 1.5e6, 3e6, 9e6])
# ts_mean  = tau_auto *0.
# dirname  = ['1_sh9e2/day4122h00.',
#     '2_sh9e3/day5496h00.',
#     '3_sh9e4/day2748h00.',
#     '4_sh9e5/day5496h00.',
#     '5_sh1.5e6/day8244h00.',
#     '6_sh3e6/day17175h00.',
#     '7_sh9e6/day12366h00.']
# for i in range(len(dirname)):
#     filename = '../data/SH/'+dirname[i]+'atmos_avg.nc'
#     # print filename
#     f = netcdf.netcdf_file(filename, 'r')
#     lat = f.variables['grid_yt'].data
#     lon = f.variables['grid_xt'].data
#     ts  = f.variables['bucket_depth'].data[0,:,:]
#     ts_mean[i] = area_mean(ts, lat)
#     f.close()
# fig = plt.figure(figsize=(4,3))
# cmap = mpl.cm.RdBu_r 
# ax = fig.add_subplot(111)
# ax.semilogx(tau_auto/86400, ts_mean,'k-',marker='s',label='SP')
# print (ts_mean)

# # tau_auto = np.array([9e2, 9e3, 9e4, 9e5,1.5e6, 3e6, 9e6])
# dirname  = ['1_eq9e2/day5496h00.',
#             '2_eq9e3/day5496h00.',
#             '3_eq9e4/day2748h00.',
#             '4_eq9e5/day5496h00.',
#             '5_eq1.5e6/day10992h00.',
#             '6_eq3e6/day10305h00.',
#             '7_eq9e6/day7557h00.']
# ts_mean  = tau_auto *0.
# for i in range(len(dirname)):
#     filename = '../data/EQ/'+dirname[i]+'atmos_avg.nc'   
#     f = netcdf.netcdf_file(filename, 'r')
#     lat = f.variables['grid_yt'].data
#     lon = f.variables['grid_xt'].data
#     ts  = f.variables['bucket_depth'].data[0,:,:]
#     ts_mean[i] = area_mean(ts, lat)
#     f.close()
# ax.semilogx(tau_auto/86400, ts_mean,'r',linestyle='--',marker='o',label='EQ')
# ax.set_xlabel(r'Conversion timescale $\tau_c$ (day)')
# ax.set_ylabel(r'Global mean T$_s$ (K)')
# ax.grid()
# ax.legend(loc=0,numpoints=1,fontsize=11)
# ax.set_xlim([1e-2,2e2])
# ax.tick_params(which='both',direction='out')
# print(ts_mean)

# plt.tight_layout()
# fig.savefig('TS_TAU_EQ.pdf')
filename = '../data/SH/1_sh9e2/day4122h00.atmos_static.nc'
f = netcdf.netcdf_file(filename, 'r')
lat = f.variables['grid_yt'].data
lon = f.variables['grid_xt'].data
zsurf  = f.variables['zsurf'].data[:,:]
f.close()

water = zsurf*0.
fig = plt.figure(figsize=(4,5))
cmap = mpl.cm.RdBu_r 
ax = fig.add_subplot(211)
X,Y = np.meshgrid(lon, lat)

water[Y<-80] = 1.
CS2 = ax.contour(X,Y,zsurf/1e3,colors='k',linewidths=0.4)
ax.contourf(X,Y,water,levels=[0,0.5,1],colors=['sandybrown','w'])#,zorder=2.5)
ax.tick_params(which='both',direction='out')
ax.set_ylabel('Latitude (deg)')
# ax.text(300, 50, 'bare ground')
ax.text(0.02, 1-0.03, 'SP', fontsize=14, 
    ha="left", va="top",
    bbox=dict(boxstyle="square",
        ec='k',
        fc='w'),
    transform=ax.transAxes)
 
water = zsurf*0.
ax = fig.add_subplot(212)
for i in range(len(lat)):
    for j in range(len(lon)):
        if (zsurf[i,j]>4e3 and abs(lat[i])<25):
            water[i,j] = 1.
CS2 = ax.contour(X,Y,zsurf/1e3,colors='k',linewidths=0.4)
ax.contourf(X,Y,water,levels=[0,0.6,1],colors=['sandybrown','w'])#,zorder=2.5)
ax.tick_params(which='both',direction='out')
ax.set_ylabel('Latitude (deg)')
ax.set_xlabel('Longitude (deg)')
# ax.text(0, 1.05, '(b) EQ', fontsize=12, transform=ax.transAxes)
ax.text(0.02, 1-0.03, 'EQ', fontsize=14, 
    ha="left", va="top",
    bbox=dict(boxstyle="square",
        ec='k',
        fc='w'),
    transform=ax.transAxes)
# CB = fig.colorbar(CS2, ax=ax,
#     orientation='horizontal')#es.ravel().tolist(), pad=0.1, location='bottom')
# CB.ax.set_xlabel(r'Water cloud mass concentration (kg/kg)')

plt.show()
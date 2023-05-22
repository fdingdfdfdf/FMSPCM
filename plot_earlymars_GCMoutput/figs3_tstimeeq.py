# -*- coding: utf-8 -*-

from turtle import title
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
mass_to_depth = 3./4*1e2 #from cloud mass to optical thickness
fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True)
cmap = mpl.cm.RdBu_r 
titles = [r'(a) $\tau_c$=0.01 day',r'(b) $\tau_c$=10 day',r'(c) $\tau_c$=100 day',
    r'(d) $\tau_c$=0.01 day',r'(e) $\tau_c$=10 day',r'(f) $\tau_c$=100 day']


tau_auto = np.array([9e2, 9e3, 9e4, 9e5, 1.5e6, 3e6, 9e6])
ts_mean  = tau_auto *0.
rain_mean  = tau_auto *0.
dirname  = ['1_sh9e2/day4122h00.',
    '4_sh9e5/day5496h00.',
    # '5_sh1.5e6/day8244h00.',
    # '6_sh3e6/day17175h00.',
    '7_sh9e6/day12366h00.']
for i in range(len(dirname)):
    filename = '../data/SH/'+dirname[i]+'atmos_daily.nc'
    # print filename
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    time= f.variables['time'].data
    ts  = f.variables['ts'].data[:,4,:] #81S
    sflux  = f.variables['flux_lhe'].data[:,:,:]
    # rain   = f.variables['total_rain'].data[0,:,:]
    # rain_mean[i]  = area_mean(rain, lat)
    # ts_mean[i] = area_mean(ts, lat)
    f.close()

#plot
    ax = axes.flat[i*2]
    ax.plot(time-time[0], np.mean(ts,axis=-1),'k-',label='SP')#r'1bar, 25$^\circ$OBL')

    ax = axes.flat[i*2+1]
    sflux_t = time*1.
    for jj in range(len(time)):
        sflux_t[jj] = area_mean(sflux[jj,:,:], lat)
    ax.plot(time-time[0], sflux_t,'k-',label='SP')#r'1bar, 25$^\circ$OBL')

# ts_mean  = tau_auto *0.
# rain_mean  = tau_auto *0.
# dirname  = ['1_obl9e2/day5496h00.',
#     # '2_obl9e3/day5496h00.',
#     # '3_obl9e4/day2748h00.',
#     '4_obl9e5/day5496h00.',
#     # '5_obl1.5e6/day10992h00.',
#     # '6_obl3e6/day10992h00.',
#     '7_obl9e6/day10992h00.']
# for i in range(len(dirname)):
#     filename = '../data/OBL45SH/'+dirname[i]+'atmos_daily.nc'
#     # print filename
#     f = netcdf.netcdf_file(filename, 'r')
#     lat = f.variables['grid_yt'].data
#     lon = f.variables['grid_xt'].data
#     time= f.variables['time'].data
#     ts  = f.variables['ts'].data[:,4,:]
#     # rain   = f.variables['total_rain'].data[0,:,:]
#     # rain_mean[i]  = area_mean(rain, lat)
#     # ts_mean[i] = area_mean(ts, lat)
#     f.close()

# #plot
#     ax = axes.flat[i]
#     ax.plot(time-time[0], np.mean(ts,axis=-1),'r--',label=r'1bar, 45$^\circ$OBL')

# tau_auto = np.array([9e2, 9e3, 9e4, 9e5,1.5e6, 3e6, 9e6])
# dirname  = ['1_ps9e2/day5496h00.',
# #             '2_ps9e3/day5496h00.',
# #             '3_ps9e4/day2748h00.',
#             '4_ps9e5/day5496h00.',
# #             '5_ps1.5e6/day16488h00.',
# #             '6_ps3e6/day16488h00.',
#             '7_ps9e6/day13740h00.']
# ts_mean  = tau_auto *0.
# rain_mean  = tau_auto *0.
# for i in range(len(dirname)):
#     filename = '../data/PS0.5SH/'+dirname[i]+'atmos_daily.nc'   
#     f = netcdf.netcdf_file(filename, 'r')
#     lat = f.variables['grid_yt'].data
#     lon = f.variables['grid_xt'].data
#     time= f.variables['time'].data
#     ts  = f.variables['ts'].data[:,4,:]
#     # rain   = f.variables['total_rain'].data[0,:,:]
#     # rain_mean[i]  = area_mean(rain, lat)
#     # ts_mean[i] = area_mean(ts, lat)
#     f.close()
dirname  = ['1_eq9e2/day5496h00.',
            '4_eq9e5/day5496h00.',
            '7_eq9e6/day7557h00.']
for i in range(len(dirname)):
    filename = '../data/EQ/'+dirname[i]+'atmos_daily.nc'   
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    time= f.variables['time'].data
    ts  = f.variables['ts'].data[:,46,37]
    sflux  = f.variables['flux_lhe'].data[:,:,:]
    f.close()
#plot
    ax = axes.flat[i*2]
    ax.plot(time-time[0], ts,'b:',label='EQ')#r'0.5bar, 25$^\circ$OBL')
    ax.set_ylabel(r'T$_s$ (K)')
    ax.tick_params(which='both',direction='out')
    ax.text(0, 1.05, titles[i],transform=ax.transAxes, fontsize=12)

    ax = axes.flat[i*2+1]
    sflux_t = time*1.
    for jj in range(len(time)):
        sflux_t[jj] = area_mean(sflux[jj,:,:], lat)
    ax.plot(time-time[0], sflux_t,'b:',label='EQ')#r'1bar, 25$^\circ$OBL')
    ax.set_ylabel(r'F$_{evap}$ (W m$^{-2}$)')
    ax.tick_params(which='both',direction='out')
    ax.text(0, 1.05, titles[i+3],transform=ax.transAxes, fontsize=12)

print lat[4], lat[46]
ax.set_xlabel(r'Time (day)')
axes.flat[-2].set_xlabel(r'Time (day)')
# ax.grid()
axes.flat[0].legend(numpoints=1, loc=0, fontsize=10)
# # ax.legend(numpoints=1,frameon=False,
# #     bbox_to_anchor=(0.3, -0.5, 1.3, .102), loc='lower left',
# #     ncol=2, mode="expand", borderaxespad=0.)
# ax.set_xlim([1e-2,2e2])
# print(ts_mean)
# 
# plt.tight_layout()
# fig.savefig('TS_TAU_EQ.pdf')

plt.show()
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

tau_auto = np.array([9e2, 9e3, 9e4, 9e5, 1.5e6, 3e6, 9e6])
ts_mean  = tau_auto *0.
rain_mean  = tau_auto *0.
dirname  = ['1_sh9e2/day4122h00.',
    '2_sh9e3/day5496h00.',
    '3_sh9e4/day2748h00.',
    '4_sh9e5/day5496h00.',
    '5_sh1.5e6/day8244h00.',
    '6_sh3e6/day17175h00.',
    '7_sh9e6/day12366h00.']
for i in range(len(dirname)):
    filename = '../data/SH/'+dirname[i]+'atmos_avg.nc'
    # print filename
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    ts  = f.variables['ts'].data[0,:,:]
    rain   = f.variables['total_rain'].data[0,:,:]
    rain_mean[i]  = area_mean(rain, lat)
    ts_mean[i] = area_mean(ts, lat)
    f.close()
fig = plt.figure(figsize=(8,4))
cmap = mpl.cm.RdBu_r 

# ax = fig.add_subplot(223)
# ax2= fig.add_subplot(224) 
# ax.semilogx(tau_auto/86400, ts_mean,'k-',marker='s',label=r'r=10 $\mu$m')
# print (rain_mean)
# ax2.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
#     'k-',marker='s',label=r'r=10 $\mu$m')

ax = fig.add_subplot(121)
ax2= fig.add_subplot(122) 
# ax.semilogx(tau_auto/86400, ts_mean,'k-',marker='s',label=r'SP')
# # print (ts_mean)
# ax2.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
#     'k-',marker='s',label=r'SP')

#KSMR21, warm start
tau_auto = np.array([2e6, 4e6, 9e6, 2e7])
dirname  = ['1_kite2e6/day6870h00.',
           '3_kite4e6warm/day13053h00.',
           '5_kite9e6warm/day10992h00.',
           '7_kite2e7warm/day19236h00.',]
ts_mean  = tau_auto *0.
rain_mean  = tau_auto *0.
for i in range(len(dirname)):
    filename = '../data/KITE21/'+dirname[i]+'atmos_avg.nc'   
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    ts  = f.variables['ts'].data[0,:,:]
    rain   = f.variables['total_rain'].data[0,:,:]
    rain_mean[i]  = area_mean(rain, lat)
    ts_mean[i] = area_mean(ts, lat)
    f.close()

ax.semilogx(tau_auto/86400, ts_mean,'r',linestyle='-',marker='o',label=r'KSMR21, warm start')
ax2.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
    'r-',marker='o',label=r'KSMR21, warm start')


#KSMR21, more ice
tau_auto = np.array([4e6])
dirname  = ['8_kite4e6iceplus/day2748h00.']
ts_mean  = tau_auto *0.
rain_mean  = tau_auto *0.
for i in range(len(dirname)):
    filename = '../data/KITE21/'+dirname[i]+'atmos_avg.nc'   
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    ts  = f.variables['ts'].data[0,:,:]
    rain   = f.variables['total_rain'].data[0,:,:]
    rain_mean[i]  = area_mean(rain, lat)
    ts_mean[i] = area_mean(ts, lat)
    f.close()

ax.semilogx(tau_auto/86400, ts_mean,'k',linestyle='none',marker='s',label=r'KSMR21, more surface ice')
ax2.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
    'k-',marker='s',label=r'KSMR21, more surface ice')

#KSMR21, cold start
tau_auto = np.array([2e6, 4e6, 9e6, 2e7])
dirname  = ['1_kite2e6/day6870h00.',
           '2_kite4e6cold/day4809h00.',
           '4_kite9e6cold/day7557h00.',
           '6_kite2e7cold/day13053h00.',]
ts_mean  = tau_auto *0.
rain_mean  = tau_auto *0.
for i in range(len(dirname)):
    filename = '../data/KITE21/'+dirname[i]+'atmos_avg.nc'   
    f = netcdf.netcdf_file(filename, 'r')
    lat = f.variables['grid_yt'].data
    lon = f.variables['grid_xt'].data
    ts  = f.variables['ts'].data[0,:,:]
    rain   = f.variables['total_rain'].data[0,:,:]
    rain_mean[i]  = area_mean(rain, lat)
    ts_mean[i] = area_mean(ts, lat)
    f.close()
ax.semilogx(tau_auto/86400, ts_mean,'b',linestyle='-',marker='^',label=r'KSMR21, cold start')
ax.set_xlabel(r'Cloud lifetime $\tau_c$ (day)')
ax.set_ylabel(r'Global mean T$_s$ (K)')
ax.grid()
# ax.legend(numpoints=1, loc='center right', fontsize=10)
ax.legend(numpoints=1,frameon=True,
    bbox_to_anchor=(0., -0.6, 2.3, .102), loc='lower left',
    ncol=2, mode="expand", borderaxespad=0.)
ax.set_xlim([10,500])
ax.tick_params(which='both',direction='out')
ax.text(0, 1.05, '(a) Surface temperature',transform=ax.transAxes, fontsize=12)

ax2.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
    'b-',marker='^',label=r'KSMR21, cold start')
ax2.plot([10, 5e2], [1,1], color='y', linewidth=5, alpha=0.5)
ax2.set_xlabel(r'Cloud lifetime $\tau_c$ (day)')
ax2.set_ylabel('Global mean \n cloud optical thickness')
ax2.grid()
# ax2.legend(loc='upper left',numpoints=1, fontsize=10)
ax2.set_xlim([10,5e2])
ax2.tick_params(which='both',direction='out')
ax2.text(0, 1.05, '(b) Cloud optical thickness (COT)',transform=ax2.transAxes, fontsize=12)
ax2.text(0.8, 0.65, 'COT=1', transform=ax2.transAxes,
    color='y',fontsize=10)

# # plt.tight_layout()
# # fig.savefig('TS_TAU_EQ.pdf')
# mass_to_depth = 3./4 /(1e3*20e-6) #from cloud mass to optical thickness
# ts_mean  = tau_auto *0.
# rain_mean  = tau_auto *0.
# dirname  = ['1_re9e2/day6870h00.',
#     '2_re9e3/day6870h00.',
#     '3_re9e4/day4122h00.',
#     '4_re9e5/day6183h00.',
#     '5_re1.5e6/day6870h00.',
#     '6_re3e6/day6870h00.',
#     '7_re9e6/day20610h00.']
# for i in range(len(dirname)):
#     filename = '../data/RE20SH/'+dirname[i]+'atmos_avg.nc'
#     # print filename
#     f = netcdf.netcdf_file(filename, 'r')
#     lat = f.variables['grid_yt'].data
#     lon = f.variables['grid_xt'].data
#     ts  = f.variables['ts'].data[0,:,:]
#     rain   = f.variables['total_rain'].data[0,:,:]
#     rain_mean[i]  = area_mean(rain, lat)
#     ts_mean[i] = area_mean(ts, lat)
#     f.close()

# ax = fig.add_subplot(223)
# ax2= fig.add_subplot(224) 
# ax.semilogx(tau_auto/86400, ts_mean,'r--',marker='o',label=r'r=20 $\mu$m')
# print (rain_mean)
# ax2.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
#     'r--',marker='o',label=r'r=20 $\mu$m')
# ax.set_xlabel(r'Cloud lifetime $\tau_c$ (day)')
# ax.set_ylabel(r'Global mean T$_s$ (K)')
# ax.grid()
# ax.legend(numpoints=1, loc='upper left', fontsize=10)
# # ax.legend(numpoints=1,frameon=False,
# #     bbox_to_anchor=(0.3, -0.5, 1.3, .102), loc='lower left',
# #     ncol=2, mode="expand", borderaxespad=0.)
# ax.set_xlim([1e-2,2e2])
# ax.tick_params(which='both',direction='out')
# ax.text(0, 1.05, '(c) Surface temperature',transform=ax.transAxes, fontsize=12)

# ax2.plot([1e-2, 2e2], [1,1], color='y', linewidth=5, alpha=0.5)
# ax2.set_xlabel(r'Cloud lifetime $\tau_c$ (day)')
# ax2.set_ylabel('Global mean \n cloud optical thickness')
# ax2.grid()
# ax2.legend(loc='upper left',numpoints=1, fontsize=10)
# ax2.set_xlim([1e-2,2e2])
# ax2.tick_params(which='both',direction='out')
# ax2.text(0, 1.05, '(d) Cloud optical thickness (COT)',transform=ax2.transAxes, fontsize=12)
# ax2.text(0.5, 0.62, 'COT=1', transform=ax2.transAxes,
#     color='y',fontsize=10)

plt.show()
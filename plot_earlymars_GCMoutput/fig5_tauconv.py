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
rco2 = 8.31/44e-3 
rho_ice = 0.917e3

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

def tauconv(slope):
    x1 = slope * 200 #slope in um
    t1 = 3*60 # in sec
    return 6*t1 / (np.exp(-x1)*(x1**3+3* x1**2+6* x1+6))

def tsed(temp, pres, radius):
    eta = 10e-6
    rho = pres/rco2/ temp
    vm  = (8*rco2*temp/np.pi)**0.5
    mfp = 2*eta / (rho*vm)
    kn  = mfp / (radius*1e-6)
    al  = 1.246+0.42*np.exp(-0.87/kn)

    vt = 2./9 * rho_ice * grav / eta *(radius*1e-6)**2 *(1+al*kn)
    ht = rco2 *temp/grav
    return ht/vt

re = np.linspace(5, 40, 100)
slope = 3./2/re
t1 = tauconv(slope)

vt = 16e3*0.5 *grav /(3*9.06**2*0.292*17e-6) *(re*1e-6)**2
ht = rco2 * 250./grav
t2 = ht/ vt

temp1 = 170.; pres1 = 0.2e5
temp2 = 160.; pres2 = 0.1e5
t2 = tsed(temp1, pres1, re)

# # print 16e3*9.8 /(3*9.06**2*0.292*17e-6) /1e8
# tau_auto = np.array([9e2, 9e5, 9e6])
# sflux_mean  = tau_auto *0.
# rain_mean   = tau_auto *0.
# cld_mean    = tau_auto *0.
# cond_mean = tau_auto*0.; cevap_mean = tau_auto*0.
# dirname  = ['1_sh9e2/day4122h00.',
#     # '2_sh9e3/day5496h00.',
#     # '3_sh9e4/day2748h00.',
#     '4_sh9e5/day5496h00.',
#     # '5_sh1.5e6/day8244h00.',
#     # '6_sh3e6/day17175h00.',
#     '7_sh9e6/day12366h00.']
# for i in range(len(dirname)):
#     filename = '../data/SH/'+dirname[i]+'atmos_avg.nc'
#     f = netcdf.netcdf_file(filename, 'r')
#     lat = f.variables['grid_yt'].data
#     lon = f.variables['grid_xt'].data
#     sflux  = f.variables['flux_lhe'].data[0,:,:]
#     rain   = f.variables['total_rain'].data[0,:,:]
#     cloud = f.variables['cld_amt'].data[0,:,:,:] #pres, lat, lon
#     cond  = f.variables['total_cond'].data[0,:,:,:] #pres, lat, lon
#     cevap = f.variables['cld_evap'].data[0,:,:,:] #pres, lat, lon

#     sflux_mean[i] = area_mean(sflux, lat)
#     rain_mean[i]  = area_mean(rain, lat)

#     ps    = f.variables['ps'].data[0,:,:]
#     phalf = np.zeros((len(pk),len(lat),len(lon)))
#     pres  = cloud *1.
#     for ij in range(len(pk)):
#         phalf[ij,:,:] = ps[:,:]*bk[ij] + pk[ij]
#     pres[:,:,:] = (phalf[1:,:,:] - phalf[:-1,:,:])  \
#         / (np.log(phalf[1:,:,:]) - np.log(phalf[:-1,:,:]))
#     cld_mean[i] = area_mean(vertical_integration(cloud, phalf), lat)
#     cond_mean[i] = area_mean(vertical_integration(cond, phalf), lat)
#     cevap_mean[i] = area_mean(vertical_integration(cevap, phalf), lat)

#     f.close()
# mass_to_depth = 3./4*1e2 #from cloud mass to optical thickness

# rain_mean = rain_mean /Lv
# fig = plt.figure(figsize=(8,3))
# cmap = mpl.cm.RdBu_r 
# titles = [r'SP, $\tau_c$=0.01 day',r'SP, $\tau_c$=10 day',r'SP, $\tau_c$=100 day']
# labels = ['Cond','Precip','Evap']
# ax = fig.add_subplot(131)
# ax.pie(np.array([cond_mean[0], rain_mean[0], cevap_mean[0]])/cond_mean[0] *100, 
#     labels=labels, autopct='%.0f%%', shadow=False)
# print cond_mean[0], cevap_mean[0], rain_mean[0]
# ax.text(0, 1.05, '(a) '+titles[0], transform=ax.transAxes, fontsize=12)
# ax.text(0,-0.15, r'Cond=2.19e-9 kg m$^{-2}$ s$^{-1}$',
#     transform=ax.transAxes, fontsize=12)

# ax = fig.add_subplot(132)
# ax.pie(np.array([cond_mean[1], rain_mean[1], cevap_mean[1]])/cond_mean[1] *100, 
#     labels=labels, autopct='%.0f%%', shadow=False)
# print cond_mean[1], cevap_mean[1], rain_mean[1]
# ax.text(0, 1.05, '(b) '+titles[1], transform=ax.transAxes, fontsize=12)
# ax.text(0,-0.15, r'Cond=1.16e-8 kg m$^{-2}$ s$^{-1}$',
#     transform=ax.transAxes, fontsize=12)

# ax = fig.add_subplot(133)
# ax.pie(np.array([cond_mean[2], rain_mean[2], cevap_mean[2]])/cond_mean[2] *100, 
#     labels=labels, autopct='%.0f%%', shadow=False)
# print cond_mean[2], cevap_mean[2], rain_mean[2]
# ax.text(0, 1.05, '(c) '+titles[2], transform=ax.transAxes, fontsize=12)
# ax.text(0,-0.15, r'Cond=1.29e-6 kg m$^{-2}$ s$^{-1}$',
#     transform=ax.transAxes, fontsize=12)

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
# ax.semilogy(re, (t1*t2)/(t1+t2) / 86400, 'k-',linewidth = 2, label='Conversion timescale')
# ax.semilogy(re, t1/ 86400, 'r--', label='Cloud ice to snow')
ax.semilogy(re, t2/ 86400, 'k-', label='Gravitational settling')
ax.plot([10, 10], [1.5e6/86400,3e6/86400], color='y', linewidth=10, alpha=0.5)
ax.plot([20, 20], [3e6/86400,9e6/86400], color='y', linewidth=10, alpha=0.5)
# ax.semilogy(10, 1.5e6/86400, linestyle='none', marker='*', ms=12, 
#     label='Transition timescale')
# ax.semilogy(20, 3e6/86400, linestyle='none', marker='*', ms=12, mfc='blue')
ax.set_xlabel(r'Cloud radius ($\mu$m)')
#ax.set_ylabel(r'Global mean cloud path (kg m$^{-2}$)')
ax.set_ylabel('Timescale (day)')
ax.grid()
ax.set_ylim([0.5, 200])
# ax.legend(numpoints=1,frameon=False,
#     bbox_to_anchor=(-0.3, -0.5, 1.5, .102), loc='lower left',
#     ncol=2, mode="expand", borderaxespad=0.)
ax.legend(numpoints=1, loc='lower left', fontsize=10)
# ax.set_ylim([1e-2,1e4])
ax.tick_params(which='both',direction='out')
# ax.text(0, 1.05, '(d)',transform=ax.transAxes, fontsize=12)
# print sflux_mean, rain_mean
# print area_mean(ps, lat[:])
# q_mean = cld_mean / 

# ax.loglog(tau_auto/86400, sflux_mean/Lv*tau_auto* mass_to_depth,
#     'r',linestyle='none',marker='o',label='Evaporation, EQ')
# # ax.loglog(tau_auto/86400, rain_mean/Lv*tau_auto* mass_to_depth, 
# #     '--',label='Precipitation')
# ax.loglog(tau_auto/86400, cld_mean* mass_to_depth, 'r--', 
#     label='Cloud mass, EQ')
# ax.plot([1e0, 2e2], [1,1], color='y', linewidth=5, alpha=0.5)
# ax.set_xlabel(r'Conversion timescale $\tau_c$ (day)')
# #ax.set_ylabel(r'Global mean cloud path (kg m$^{-2}$)')
# ax.set_ylabel('Global mean \n cloud optical thickness')
# ax.grid()
# # ax.legend(numpoints=1,frameon=False,
# #     bbox_to_anchor=(-0.3, -0.5, 1.5, .102), loc='lower left',
# #     ncol=2, mode="expand", borderaxespad=0.)
# ax.legend(numpoints=1, loc='upper left', fontsize=10)
# ax.set_xlim([1e-2,2e2])
# ax.tick_params(which='both',direction='out')
# ax.text(0, 1.05, '(a) GCM',transform=ax.transAxes, fontsize=12)

# print sflux_mean, rain_mean

# # plt.tight_layout()
# # fig.savefig('SFLUX_CLD_EQ.pdf')


plt.show()
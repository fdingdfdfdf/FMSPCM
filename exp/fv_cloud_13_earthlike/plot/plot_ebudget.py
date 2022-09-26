# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font',family = 'serif',size = 11)

from scipy.io import netcdf_file
from scipy.interpolate import interp1d
# import netCDF4

rp  = 6371e3
rcp = 8.31/32e-3 / 1062.5
spin= 7.292e-5 
grav = 9.8
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

# # filename = '../day7602h00.atmos_avg.nc'
filename = '../earth_mean.nc'
f = netcdf_file(filename, 'r')
print(f.history)
lat = f.variables['grid_yt'].data
lon = f.variables['grid_xt'].data
# # # pfull = f.variables['pfull'].data
# # # uwind = np.mean(f.variables['ucomp'].data[0,:,:,:], axis=2)
# # # olrnu = f.variables['OLRnu'].data[0,:,:,:]
# # ts    = f.variables['ts'].data[0,:,:]
# flux_t   = f.variables['flux_t'].data[0,:,:]
# flux_lhe = f.variables['flux_lhe'].data[0,:,:]
# ps    = f.variables['ps'].data[0,:,:]
# rain  = f.variables['total_rain'].data[0,:,:]
# temp  = f.variables['temp'].data[0,:,:,:] #pres, lat, lon
# uwind = f.variables['ucomp'].data[0,:,:,:] #pres, lat, lon
# # omega = f.variables['omega'].data[0,:,:,:] #pres, lat, lon
# # vwind = f.variables['vcomp'].data[0,:,:,:] #pres, lat, lon
# # sphum = f.variables['sphum'].data[0,:,:,:] #pres, lat, lon
# cld   = f.variables['cld_amt'].data[0,:,:,:] #pres, lat, lon
f.close()
# 
# filename = '../earth_e.nc'
# f = netcdf_file(filename, 'r')
# olr = f.variables['OLR'].data[0,:,:]
# osr = f.variables['OSR'].data[0,:,:]
# isr = f.variables['solar'].data[0,:,:]
# print(f.history)
# f.close()
# 
# # #heating rate
# filename = '../day9490h00/day9490h00.atmos_lbl.nc'
# f     = netcdf_file(filename, 'r')
# # hrlw  = f.variables['hrlw'].data[0,:,:,:]
# # hrsw  = f.variables['hrsw'].data[0,:,:,:]
# # hr    = hrlw[::-1,:,:] + hrsw[::-1,:,:]
# glr = f.variables['flux_lw'].data[0,0,:,:] #upward, note the index
# gsr = f.variables['flux_sw'].data[0,0,:,:]
# f.close()

filename = '../day9490h00/day9490h00.spectrum.nc'
f     = netcdf_file(filename, 'r')
# hrlw  = f.variables['hrlw'].data[0,:,:,:]
# hrsw  = f.variables['hrsw'].data[0,:,:,:]
# hr    = hrlw[::-1,:,:] + hrsw[::-1,:,:]
olrnu = f.variables['OLRnu'].data[0,:,92/2,:] #upward, note the index
solarnu = f.variables['solarnu'].data[0,:,:,:] #wavenumber, lat, lon
f.close()

# #for vertical interpolation
# filename = '../day9490h00/day9490h00.atmos_static.tile6.nc'
# f     = netcdf_file(filename, 'r')
# pk    = f.variables['pk'].data
# bk    = f.variables['bk'].data
# f.close()
# phalf = np.zeros((len(pk),len(lat),len(lon)))
# pres  = temp *1.
# for i in range(len(pk)):
#     phalf[i,:,:] = ps[:,:]*bk[i] + pk[i]
# cld1 = np.sum(cld[:,:,:]*(phalf[1:,:,:]-phalf[:-1,:,:])/grav,  axis=0)

# pres[:,:,:] = (phalf[1:,:,:] - phalf[:-1,:,:])  \
#     / (np.log(phalf[1:,:,:]) - np.log(phalf[:-1,:,:]))

# # #compute geo height
# # gph = temp *0.; gph_h = phalf *0.
# # for i in range(len(pk)-2,-1,-1):
# #     gph_h[i,:,:] = gph_h[i+1,:,:] + ro2 * temp[i,:,:] \
# #         * (np.log(phalf[i+1,:,:]) - np.log(phalf[i,:,:]))
# #     gph[i,:,:] = gph_h[i+1,:,:] + ro2 * temp[i,:,:] \
# #         * (np.log(phalf[i+1,:,:]) - np.log(pres[i,:,:]))    

# # levelp = np.logspace(np.log10(200), 5, 30) 
# # levelp = np.logspace(2, 5, 30)  #ps=1bar
# levelp = np.logspace(2, np.log10(98000), 30)
# y1in   = temp *1.; y1out  = np.zeros((30,len(lat),len(lon)))
# y2in   = cld *1.; y2out  = np.zeros((30,len(lat),len(lon)))
# y3in   = vwind *1.; y3out  = np.zeros((30,len(lat),len(lon)))
# y4in   = sphum *1.; y4out  = np.zeros((30,len(lat),len(lon)))

# for i in range(len(lat)):
#     for j in range(len(lon)):
#         ff = interp1d(np.log(pres[:,i,j]), y1in[:,i,j], 
#              bounds_error=False, fill_value=np.nan)
#         y1out[:,i,j] = ff(np.log(levelp))
#         ff = interp1d(np.log(pres[:,i,j]), y2in[:,i,j], 
#              bounds_error=False, fill_value=np.nan)
#         y2out[:,i,j] = ff(np.log(levelp))
#         ff = interp1d(np.log(pres[:,i,j]), y3in[:,i,j], 
#              bounds_error=False, fill_value=np.nan)
#         y3out[:,i,j] = ff(np.log(levelp))
#         ff = interp1d(np.log(pres[:,i,j]), y4in[:,i,j], 
#              bounds_error=False, fill_value=np.nan)
#         y4out[:,i,j] = ff(np.log(levelp))

# slat = np.sin(lat / 180.*np.pi) 
# tlat = np.tan(lat / 180.*np.pi) 
# clat = np.cos(lat / 180.*np.pi) 

# rh = y4out *1.
# for i in range(len(levelp)):
#     rh[i,:,:] = levelp[i]*y4out[i,:,:]/(0.622+0.378*y4out[i,:,:]) \
#         / esat(y1out[i,:,:])
# rh_m = np.mean(rh,axis=2)

# #stream function
# vm = np.mean(y3out, axis=2)
# stf = np.zeros(vm.shape)
# for i in range(len(levelp)-1):
#     if (any(np.isnan(vm[i,:]))):
#         continue
#     else:
#         stf[i+1,:] = stf[i,:] + (vm[i,:]+vm[i+1,:])/2. \
#         * (levelp[i+1] - levelp[i])/g *2*np.pi*rp*clat[:]


# # ps_night = area_mean(ps[:, 0:len(lon)/2], lat[:])
# # ps_day   = area_mean(ps[:, len(lon)/2:len(lon)], lat[:])

# # plot 
# X,Y = np.meshgrid(lat, levelp)
# # # 
# #colors = ['salmon','darkorange','darkgreen','royalblue']
# # levels = range(240, 640, 15)
# # norm = mpl.cm.colors.Normalize(vmax=np.nanmax(abs(Z)), 
# #     vmin= - np.nanmax(abs(Z)))
# cmap = mpl.cm.bwr #RdBu_r

# fig = plt.figure(figsize=(8,3))
# #
# # levels = range(-130, 200, 10)
# # norm   = mpl.cm.colors.Normalize(-190, 190)
# # norm   = mpl.cm.colors.Normalize(vmax=np.nanmax(abs(Z)), 
# #     vmin= - np.nanmax(abs(Z)))
# cmap = mpl.cm.bwr #RdBu_r 
# ax = fig.add_subplot(121)
# levels = np.logspace(-15,-4,10)
# norm   = mpl.cm.colors.LogNorm(1e-15,1e-4)
# CS2 = ax.contourf(X,Y, np.mean(y2out,axis=2),
#     levels=levels, norm=norm,
#     cmap=cmap  )
# # CS2 = ax.contourf(X,Y, rh_m, levels=np.linspace(0,1,10),
# #     cmap=cmap  )
# # CS = ax.contour(X,Y,np.mean(y2out,axis=2), 
# #     levels, colors='w', linewidths=0.5)
# # CS2= ax.contourf(X,Y,np.mean(y2out,axis=2), levels,
# #     cmap=cmap, norm=norm   )
# CB = fig.colorbar(CS2, orientation='horizontal')
# CB.ax.set_xlabel(r'Cloud water (kg/kg)')

# ax.set_xlabel('Lat (deg)')
# ax.set_ylabel(r'Pressure (Pa)')
# ax.set_yscale('log')
# ax.invert_yaxis()
# ax.grid()
# ax.set_ylim([1e5, 1e2])
# # 

# # levels = np.zeros(25)
# # levels[13:25] = np.logspace(8,12,12)
# # levels[:12] = -levels[24:12:-1]

# ax = fig.add_subplot(122)
# levels = np.logspace(-6,1,10)
# norm   = mpl.cm.colors.LogNorm(1e-6,10)
# CS2 = ax.contourf(X,Y, np.mean(y2out /y4out,axis=2),
#     levels=levels, norm=norm,
#     cmap=cmap  )
# CB = fig.colorbar(CS2, orientation='horizontal')
# CB.ax.set_xlabel(r'Cloud water /water vapor(kg/kg)')
# ax.set_xlabel('Lat (deg)')
# ax.set_ylabel(r'Pressure (Pa)')
# ax.set_yscale('log')
# ax.invert_yaxis()
# ax.grid()
# ax.set_ylim([1e5, 1e2])

# # CS2 = ax.pcolor(X,Y, stf, 
# #     norm=mpl.colors.SymLogNorm(linthresh=1e7, linscale=0.01,
# #     vmin=-1e12, vmax=1e12),
# #     cmap=cmap  )
# # # CS = ax.contour(X,Y,np.mean(y2out,axis=2), 
# # #     levels, colors='w', linewidths=0.5)
# # # CS2= ax.contourf(X,Y,np.mean(y2out,axis=2), levels,
# # #     cmap=cmap, norm=norm   )
# # CB = fig.colorbar(CS2, orientation='horizontal')
# # CB.ax.set_xlabel(r'Stream function (kg/s)')
# # ax.set_xlabel('Lat (deg)')
# # ax.set_ylabel(r'Pressure (Pa)')
# # ax.set_yscale('log')
# # ax.invert_yaxis()
# # ax.grid()
# # ax.set_ylim([1e5, 1e2])

# isr_mean = area_mean(isr, lat)
# olr_mean = area_mean(olr, lat)
# print isr_mean, olr_mean
# rain_m = area_mean(rain/Lv *86400, lat)
# pgpcp_m = area_mean(p_gpcp, lat_gpcp)
# print rain_m, pgpcp_m
# 
# gsr_mean = area_mean(gsr, lat)
# glr_mean = area_mean(glr, lat)
# sh_mean  = area_mean(flux_t, lat)
# lh_mean  = area_mean(flux_lhe, lat)
# print gsr_mean, glr_mean, sh_mean, lh_mean
clat = np.cos(lat / 180.*np.pi)
var1 = np.mean(solarnu, axis=-1) *clat 
solarnu_mean = np.sum(var1, axis=-1) / np.nansum(clat)
nu_sw        = np.linspace(1e3, 5e4, 2000)

solar = np.sum((solarnu_mean[1:]+solarnu_mean[:-1])/2. \
    *(nu_sw[1:]-nu_sw[:-1]))
idx = np.abs(nu_sw - 3.33e4).argmin()
solar_uv = np.sum((solarnu_mean[idx+1:]+solarnu_mean[idx:-1])/2. \
    *(nu_sw[idx+1:]-nu_sw[idx:-1]))
print solar, solar_uv

fig = plt.figure(figsize=(5,6))
ax = fig.add_subplot(211)
labels = ['GSR', 'GLR', 'Sensible', 'Latent']
earth = [161, -63, -17, -80]
fms = [193.46, -37.01, -13.47, -144.45]
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars
rects1 = ax.bar(x - width/2, earth, width, label='Reanalysis')
rects2 = ax.bar(x + width/2, fms, width, label='FMS-PCM')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel(r'Surface energy flux (W m$^{-2}$)')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
ax.grid()
# ax.bar_label(rects1, padding=3)
# ax.bar_label(rects2, padding=3)


ax = fig.add_subplot(212)
ax.plot(nu_sw, solarnu_mean,'k-',lw=0.5,label='FMS-PCM')
ax.axvspan(nu_sw[idx], nu_sw[-1], color='red', alpha=0.5)
# ax.plot(wnum,  np.mean(olr_airs,axis=-1),'r-',lw=0.5,label='AIRS')
# ax.set_xlim([0,2500])
ax.set_xlabel('Wavenumber (cm$^{-1}$)')
ax.set_ylabel(r'Spectral ISR (W m$^{-2}$ cm$^{-1}$)')
ax.grid()
ax.legend(loc=0)
ax.text(0.7,0.3,'3.77 / 340.30',size=13,
    transform=ax.transAxes)

plt.show()
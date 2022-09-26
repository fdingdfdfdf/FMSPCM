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

filename = '../day7602h00.atmos_avg.nc'
# filename = '../default_mean.nc'
f = netcdf_file(filename, 'r')
print(f.history)
lat = f.variables['grid_yt'].data
lon = f.variables['grid_xt'].data
# # pfull = f.variables['pfull'].data
# # uwind = np.mean(f.variables['ucomp'].data[0,:,:,:], axis=2)
# # olrnu = f.variables['OLRnu'].data[0,:,:,:]
# ts    = f.variables['ts'].data[0,:,:]
ps    = f.variables['ps'].data[0,:,:]
rain  = f.variables['total_rain'].data[0,:,:]
flux_t= f.variables['flux_t'].data[0,:,:]
temp  = f.variables['temp'].data[0,:,:,:] #pres, lat, lon
uwind = f.variables['ucomp'].data[0,:,:,:] #pres, lat, lon
# omega = f.variables['omega'].data[0,:,:,:] #pres, lat, lon
# vwind = f.variables['vcomp'].data[0,:,:,:] #pres, lat, lon
# sphum = f.variables['sphum'].data[0,:,:,:] #pres, lat, lon
cld   = f.variables['cld_amt'].data[0,:,:,:] #pres, lat, lon
f.close()

filename = '../day7602h00.atmos_echeck.nc'
f = netcdf_file(filename, 'r')
olr = f.variables['OLR'].data[0,:,:]
osr = f.variables['OSR'].data[0,:,:]
print(f.history)
f.close()

filename = '../day7602h00.atmos_static.nc'
f = netcdf_file(filename, 'r')
print(f.history)
isr = f.variables['solar'].data[0,:,:]
f.close()

# filename = '../obs/CERES_EBAF-TOA_Ed4.1_Subset_CLIM01-CLIM12.nc'
# f = netcdf_file(filename, 'r') 
# lat_ceres = f.variables['lat'].data
# isr_ceres = f.variables['zsolar_clim'].data[:,:] #time,lat
# osr_ceres = f.variables['ztoa_sw_all_clim'].data[:,:] #time,lat
# olr_ceres = f.variables['ztoa_lw_all_clim'].data[:,:] #time,lat
# f.close()
filename = '../obs/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
f = netcdf_file(filename, 'r') 
lat_ceres = f.variables['lat'].data
isr_ceres = f.variables['zsolar_clim'].data[:,:] #time,lat
osr_ceres = f.variables['ztoa_sw_all_clim'].data[:,:] #time,lat
olr_ceres = f.variables['ztoa_lw_all_clim'].data[:,:] #time,lat
glr_ceres = f.variables['zsfc_net_lw_all_clim'].data[:,:] #downward
gsr_ceres = f.variables['zsfc_net_sw_all_clim'].data[:,:]
f.close()


filename = '../obs/air3.nc'
f = netcdf_file(filename, 'r') 
lat_ncep = f.variables['lat'].data
lon_ncep = f.variables['lon'].data
st_ncep  = f.variables['air'].data[:,:,:]
sair_ncep= np.mean(st_ncep,axis=0)
f.close()

filename = '../obs/uwnd3.nc'
f = netcdf_file(filename, 'r') 
st_ncep  = f.variables['uwnd'].data[:,:,:]
uwnd_ncep= np.mean(st_ncep,axis=0)
f.close()

filename = '../obs/precip.mon.ltm.nc'
f = netcdf_file(filename, 'r') 
lat_gpcp = f.variables['lat'].data
st_gpcp  = f.variables['precip'].data[:,:,:]
p_gpcp   = np.mean(st_gpcp,axis=0)
f.close()

#heating rate
filename = '../day7607h00.atmos_lbl.nc'
f     = netcdf_file(filename, 'r')
# hrlw  = f.variables['hrlw'].data[0,:,:,:]
# hrsw  = f.variables['hrsw'].data[0,:,:,:]
# hr    = hrlw[::-1,:,:] + hrsw[::-1,:,:]
glr = f.variables['flux_lw'].data[0,0,:,:] #upward, note the index
gsr = f.variables['flux_sw'].data[0,0,:,:]
f.close()

#for vertical interpolation
filename = '../day7602h00.atmos_static.tile6.nc'
f     = netcdf_file(filename, 'r')
pk    = f.variables['pk'].data
bk    = f.variables['bk'].data
f.close()
phalf = np.zeros((len(pk),len(lat),len(lon)))
pres  = temp *1.
for i in range(len(pk)):
    phalf[i,:,:] = ps[:,:]*bk[i] + pk[i]
cld1 = np.sum(cld[:,:,:]*(phalf[1:,:,:]-phalf[:-1,:,:])/grav,  axis=0)

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

isr_mean = area_mean(isr, lat)
olr_mean = area_mean(olr, lat)
print isr_mean, olr_mean
rain_m = area_mean(rain/Lv *86400, lat)
pgpcp_m = area_mean(p_gpcp, lat_gpcp)
print rain_m, pgpcp_m

fig, axes = plt.subplots(2,2,figsize=(9,6),sharex=True)
# ax = axes[0,0]
# ax.plot(lat, np.mean(isr-osr,axis=-1),'k-',lw=2,label='FMS-PCM')
# ax.plot(lat_ceres, np.mean(isr_ceres-osr_ceres,axis=0),'k--',
#     lw=2,label='CERES')

# # xx = np.mean(isr_ceres,0)
# # clat = np.cos(lat_ceres/180.*np.pi)
# # print np.sum(xx*clat)/np.sum(clat)

# # ax.set_xlabel('Lat (deg)')
# ax.set_ylabel(r'Radiative flux (W m$^{-2}$)')
# ax.grid()
# ax.legend(loc=0)
# ax.text(0,1.05,'(a) Absorbed Shortwave Radiation (ASR)',size=13,
#     transform=ax.transAxes)

# ax = axes[0,1]
# ax.plot(lat, np.mean(olr,axis=-1),'k-',lw=2,label='FMS-PCM')
# ax.plot(lat_ceres, np.mean(olr_ceres,axis=0),'k--',
#     lw=2,label='CERES')
# # ax.set_xlabel('Lat (deg)')
# ax.set_ylabel(r'Radiative flux (W m$^{-2}$)')
# ax.grid()
# ax.legend(loc=0)
# ax.text(0,1.05,'(b) Outgoing Longwave Radiation (OLR)',size=13,
#     transform=ax.transAxes)

# ax = axes[1,0]
# ax.plot(lat, np.mean(temp[-1,:,:],axis=-1),'k-',lw=2,label='FMS-PCM')
# ax.plot(lat_ncep, np.mean(sair_ncep,axis=-1)+273.15,'k--',lw=2,
#     label='NCEP/NCAR')
# # ax.set_xlabel('Lat (deg)')
# ax.set_ylabel(r'Temperature (K)')
# ax.grid()
# ax.legend(loc=0)
# ax.text(0,1.05,'(c) Surface air temperature',size=13,
#     transform=ax.transAxes)

ax = axes[1,1]
# ax.plot(lat, np.mean(uwind[-1,:,:],axis=-1),'k-',lw=2,label='FMS-PCM')
# ax.plot(lat_ncep, np.mean(uwnd_ncep,axis=-1),'k--',lw=2,
#     label='NCEP/NCAR')
ax.plot(lat, np.mean(rain,axis=-1)/Lv *86400,
    'k-',lw=2,label='FMS-PCM')
ax.plot(lat, np.mean(cld1,axis=-1)/1200 *86400,
    'r--',lw=2,label=r'Mcld/time')
ax.plot(lat_gpcp, np.mean(p_gpcp,axis=-1),'k--',lw=2,
    label='GPCP')
# ax.set_xlabel('Lat (deg)')
ax.set_ylabel(r'Precipitation (mm/day)')
ax.grid()
ax.legend(loc=0)
ax.text(0,1.05,'(d) Precipitation rate, 5.59/2.69=2.08',size=13,
    transform=ax.transAxes)

# ax = axes[2,0]
# ax.plot(lat, np.mean(-slr,axis=-1),'k-',lw=2,label='FMS-PCM')
# ax.plot(lat_ceres, np.mean(slr_ceres,axis=0),'k--',
#     lw=2,label='CERES')
# # ax.set_xlabel('Lat (deg)')
# ax.set_ylabel(r'Radiative flux (W m$^{-2}$)')
# ax.grid()
# ax.legend(loc=0)
# ax.text(0,1.05,'(e) Net downward longwave flux at surface',size=13,
#     transform=ax.transAxes)

ax = axes[0,0]
ax.plot(lat, np.mean(isr-osr-gsr,axis=-1),'k-',lw=2,label='FMS-PCM, 52.16')
ax.plot(lat_ceres, np.mean(isr_ceres-osr_ceres-gsr_ceres,axis=0),'k--',
    lw=2,label='CERES, 77.58')
# ax.set_xlabel('Lat (deg)')
ax.set_ylabel(r'Radiative flux (W m$^{-2}$)')
ax.grid()
ax.legend(loc=0)
ax.text(0,1.05,'(a) Net shortwave flux absorption',size=13,
    transform=ax.transAxes)

ax = axes[0,1]
ax.plot(lat, np.mean(olr-glr,axis=-1),'k-',lw=2,label='FMS-PCM, 226.03')
ax.plot(lat_ceres, np.mean(olr_ceres+glr_ceres,axis=0),'k--',
    lw=2,label='CERES, 186.67')
# ax.set_xlabel('Lat (deg)')
ax.set_ylabel(r'Radiative flux (W m$^{-2}$)')
ax.grid()
ax.legend(loc=0)
ax.text(0,1.05,'(b) Net longwave flux divergence',size=13,
    transform=ax.transAxes)

ax = axes[1,0]
ax.plot(lat, np.mean(gsr-glr,axis=-1),'k-',lw=2,label='FMS-PCM, 173.87-(11.45)')
ax.plot(lat_ceres, np.mean(gsr_ceres+glr_ceres,axis=0),'k--',
    lw=2,label='CERES, 109.09-(26.11)')
# ax.set_xlabel('Lat (deg)')
ax.set_ylabel(r'Radiative flux (W m$^{-2}$)')
ax.grid()
ax.legend(loc=0)
ax.text(0,1.05,'(c) Net radiative flux divergence',size=13,
    transform=ax.transAxes)

# xx = np.mean(olr_ceres+slr_ceres,0)
xx = np.mean(isr_ceres-osr_ceres-gsr_ceres,0)
clat = np.cos(lat_ceres/180.*np.pi)
print np.sum(xx*clat)/np.sum(clat), area_mean(isr-osr-gsr,lat)
print area_mean(flux_t, lat)
print 2.69187/86400*Lv

plt.show()
# -*- coding: utf-8 -*-

import math
import numpy 
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font',family = 'serif',size = 12)

from scipy.io import netcdf
from scipy.interpolate import interp1d

rp  = 3389.5e3 #6371e3
# rcp = 8.31/32e-3 / 1062.5
spin= 7.292e-5 
grav= 3.71 #9.8
# ro2 = 8.31/32e-3 
Lv = 2.5e6
Rv = 461.5
Rd = 8.31/44e-3
zvir = Rv/Rd-1.
CCN = 1e5
Pi = numpy.pi

c = 2.99e8
h = 6.63e-34
kB = 1.38e-23
def planck(nu_cm, T):
    nu = 100. *c *nu_cm
    Btemp = (2*h*nu**3/c**2)/(numpy.exp(h*nu/(kB*T))  -  1.0) #W/m2/sr/Hz
    Bnu   = 1.0e2*c*Btemp  # W/m2/sr/cm^-1
    return Bnu*numpy.pi

## -----------------------------------------------------------
## BELOW: compute derivative with respect to temperature, without pi

def dPlanckdT_nu(nu,T):
    u = h*nu/(kB*T)
    u[u>500.] = 500. #To prevent overflow
    return 2.*h**2/(kB*c**2) * (nu**4/T**2) * numpy.exp(u)/( (numpy.exp(u)-1.)**2 )

# note: assume 'n' is in cm^-1!
def dPlanckdT_n(n,T,unit="cm^-1"):
    if unit=="m^-1":
        k = 1.
    elif unit=="cm^-1":
        k = 100.
    else:
        print( "(Planck_n) Error: unit not recognized!")

    return (k*c) * dPlanckdT_nu( c*n*k ,T)

#NANG = 90

#Translated into Python from bhmie.c

#This uses one-based arrays because it descends from a
#Fortran original.  The zero index isn't used. To rationalize
#the code, it should be re-written as zero-based, but that
#needs to be done cautiously since it's easy to introduce
#mistakes..
#
#Inputs:
#         x        2*Pi*rad*REFMED/lam
#         refrel   Relative index of refraction of particle
#         nang     Number of angles for phase function
#
#If called with the default nang, it returns only the scattering eff,
#absorption eff, and asymmetry factor.  If called with nang>2, it
#returns the theta array and the phase function as well.
#

#----------------------------------------------------------------
#
#6/13/2013: Replaced numpy with numpy and fixed Float type so
#           it works with numpy
#
#Change Log:
def bhmie(x,refrel,nang = 2):
    an_old = 0. + 0.j
    bn_old = 0. + 0.j
    amu = numpy.zeros(nang+1,numpy.float)
    theta = numpy.zeros(nang+1,numpy.float)
    pi = numpy.zeros(nang+1,numpy.float)
    tau = numpy.zeros(nang+1,numpy.float)
    pi0 = numpy.zeros(nang+1,numpy.float)
    pi1 = numpy.zeros(nang+1,numpy.float)

    d = numpy.zeros(3000,numpy.complex) #**Change this to a dynamical allocation
    s1 = numpy.zeros(2*nang,numpy.complex)
    s2 = numpy.zeros(2*nang,numpy.complex)
    dx = x
    y = x*refrel
    xstop = x + 4.*x**(1./3.) + 2.
    nstop = int(xstop)
    ymod = abs(y)
    
    if xstop > ymod:
        nmx = int(xstop + 15)
    else:
        nmx = int(ymod + 15)
    dang = (Pi/2.)/(nang - 1.)
    for j in range(1,nang+1):
        theta[j] = (j - 1.)*dang
        amu[j] = math.cos(theta[j])

    d[nmx] = 0. + 0.j
    nn = nmx - 1
    for n in range(1,nn+1):
        rn = nmx - n + 1
        d[nmx-n] =rn/y -1./(d[nmx-n+1] + rn/y)
        
    for j in range(1,nang+1):
        pi0[j] = 0.0
        pi1[j] = 1.0
    nn = 2*nang - 1

    for j in range(1,nn+1):
        s1[j] = complex(0.0,0.0)
        s2[j] = complex(0.0,0.0)

    psi0 = math.cos(dx)
    psi1 = math.sin(dx)
    chi0 = -math.sin(x)
    chi1 = math.cos(x)
    apsi0 = psi0
    apsi1 = psi1
    xi0 = complex(apsi0,-chi0)
    xi1 = complex(apsi1,-chi1)
    qsca = 0.0
    g = 0.0
    n = 1
#--------------------------------------
    while n - 1 - nstop < 0:
        dn = float(n)
        rn = float(n)
        fn = (2.*rn + 1.)/(rn*(rn + 1.))
        psi = (2.*dn - 1.)*psi1/dx - psi0
        apsi = psi
        chi = (2.*rn - 1.)*chi1/x - chi0
        xi = complex(apsi,-chi)
#-----------------------------------------        
        an =apsi*(d[n]/refrel+rn/x) - apsi1
        an = an/((d[n]/refrel+rn/x)*xi-xi1)
        bn = apsi*(refrel*d[n]+rn/x) - apsi1
        bn = bn/((refrel*d[n]+rn/x)*xi-xi1)
        qsca += (2*rn + 1.)*(abs(an)**2 + abs(bn)**2)
        if rn > 1:
            g += ((rn - 1.)*(rn + 1.)/rn)*(an_old*an.conjugate()+bn_old*bn.conjugate()).real +\
                 ((2.*(rn - 1.) + 1.)/((rn - 1.)*rn))*(an_old*bn_old.conjugate()).real                
        an_old = an
        bn_old = bn
#------------------------------------------
        for  j in range(1,nang+1):
            jj = 2*nang - j
            pi[j] = pi1[j]
            tau[j] = rn*amu[j]*pi[j] - (rn + 1)*pi0[j]
            p = (-1)**(n-1)
            s1[j] = s1[j]+fn*(pi[j]*an +tau[j]*bn)
            t = (-1)**n
            s2[j] = s2[j]+fn*(tau[j]*an+pi[j]*bn)
##      if(j == jj) continue; 
            if not (j == jj):
                s1[jj] = s1[jj] + fn*(pi[j]*p*an+tau[j]*t*bn)
                s2[jj] = s2[jj]+  fn*(tau[j]*t*an+pi[j]*p*bn)
  
        psi0 = psi1
        psi1 = psi
        apsi1 = psi1
        chi0 = chi1
        chi1 = chi
        xi1 = complex(apsi1,-chi1)
        n = n + 1
        rn = float(n)

        for j in range(1,nang+1):
            pi1[j] = ((2.*rn - 1.)/(rn - 1.))*amu[j]*pi[j]
            pi1[j] = pi1[j] - rn*pi0[j]/(rn - 1.)
            pi0[j] = pi[j]
#  while(n - 1 - nstop < 0);
#-------------------------------
#Returns
    qsca *= 2./x**2
    qext = (4./x**2)*s1[1].real 
    qback = (4./x**2)*abs(s1[2*nang - 1])**2
    g *= 4./(x**2*qsca)
    qabs = qext - qsca
    #
    #Compute the phase function and normalize it
    P = numpy.absolute(s1)**2 + numpy.absolute(s2)**2
    P = P[1:] #Convert it to a zero based array 
    thetaAll = numpy.array([j*dang for j in range(len(P))])
    sinthetaAll = numpy.sin(thetaAll)
    norm = sum(sinthetaAll*P*dang)
    P = 2.*P/norm #Normalize such that int P dOmega = 4Pi
    if nang > 2:
        return qabs,qsca,g,thetaAll,P
    else:
        return qabs,qsca,g

#----Utility programs for plotting curves
def MieParams(xList,refrel):
    qscaL = []
    qabsL = []
    gL = []
    for x in xList:
        qabs,qsca,g = bhmie(x,refrel)
        qabsL.append(qabs)
        qscaL.append(qsca)
        gL.append(g)
    return qabsL,qscaL,gL
    
# def esat(T):
#     return 610.78 * np.exp(-Lv/Rv*(1./T - 1./273.16))
    
# def area_mean(var2d, lat): 
#     clat = np.cos(lat / 180.*np.pi)
#     var1 = np.mean(var2d, axis=1) *clat 
#     var2 = np.sum(var1, axis=0) / np.nansum(clat)
#     return var2

# def vertical_interp(var3d, p3d, p1d): #pres, lat, lon, for visualization
#     aa, nlat, nlon = var3d.shape
#     varout = np.zeros((len(p1d),nlat,nlon))
#     for i in range(nlat):
#         for j in range(nlon): 
#             ff = interp1d(np.log(p3d[:,i,j]), var3d[:,i,j], 
#                     bounds_error=False, fill_value=np.nan)
#             varout[:,i,j] = ff(np.log(p1d))
#     return varout

# def vertical_integration(var3d, phalf3d):
#     npfull, nlat, nlon = var3d.shape
#     varout = np.zeros((nlat, nlon))
#     for i in range(npfull):
#         varout[:,:] = varout[:,:] + var3d[i,:,:]* \
#             (phalf3d[i+1,:,:] - phalf3d[i,:,:])/grav
#     return varout

r = 1.
waterice = numpy.loadtxt('waterice_refrective.txt') #wl in um, n, k
xList = 2.*Pi *r/waterice[:,0]
refRe = waterice[:,1]
refIm = waterice[:,2]
absList = xList*1.
scaList = xList*1.
gList   = xList*1.

for i in range(len(xList)):
    refrel = complex(refRe[i],refIm[i])
    qabs,qsca,g = bhmie(xList[i],refrel)
    absList[i] = qabs
    scaList[i] = qsca
    gList[i]   = g

extList = absList+scaList


# fig = plt.figure(figsize=(4,3))
# cmap = mpl.cm.RdBu_r 
# ax = fig.add_subplot(111)
# xx = numpy.linspace(1,3000) #in cm-1
# ax.plot(xx, planck(xx ,300))
wn = 1e4 / waterice[:,0] #in cm-1
planck_list = planck(wn, 200)
xx = (planck_list[:-1]+planck_list[1:])/2.*(wn[:-1] - wn[1:])
planck_sum = numpy.sum(xx)
print planck_sum, 5.67e-8*200**4
xx = (planck_list[:-1]+planck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (absList[:-1]+absList[1:])/2.
absmean_list = [numpy.sum(xx) /planck_sum] #Planck mean 

dplanck_list = dPlanckdT_n(wn,200)
xx = (dplanck_list[:-1]+dplanck_list[1:])/2.*(wn[:-1] - wn[1:])
dplanck_sum = numpy.sum(xx)
print dplanck_sum, 5.67e-8*200**3 *4 /Pi
xx = (dplanck_list[:-1]+dplanck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (1./absList[:-1]+1./absList[1:])/2.
rosselandabs_list = [dplanck_sum / numpy.sum(xx)]


fig = plt.figure(figsize=(4,3))
cmap = mpl.cm.RdBu_r 
ax = fig.add_subplot(111)
ax.semilogx(waterice[:,0], absList, 'g:', label=r'r=1$\mu$m')


r = 5.
xList = 2.*Pi *r/waterice[:,0]
for i in range(len(xList)):
    refrel = complex(refRe[i],refIm[i])
    qabs,qsca,g = bhmie(xList[i],refrel)
    absList[i] = qabs
    scaList[i] = qsca
    gList[i]   = g
ax.semilogx(waterice[:,0], absList, 'b-.', label=r'r=5$\mu$m')
xx = (planck_list[:-1]+planck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (absList[:-1]+absList[1:])/2.
absmean_list.append(numpy.sum(xx) /planck_sum)
xx = (dplanck_list[:-1]+dplanck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (1./absList[:-1]+1./absList[1:])/2.
rosselandabs_list.append(dplanck_sum / numpy.sum(xx))

r = 10.
xList = 2.*Pi *r/waterice[:,0]
for i in range(len(xList)):
    refrel = complex(refRe[i],refIm[i])
    qabs,qsca,g = bhmie(xList[i],refrel)
    absList[i] = qabs
    scaList[i] = qsca
    gList[i]   = g
ax.semilogx(waterice[:,0], absList, 'r--', label=r'r=10$\mu$m')
xx = (planck_list[:-1]+planck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (absList[:-1]+absList[1:])/2.
absmean_list.append(numpy.sum(xx) /planck_sum)
xx = (dplanck_list[:-1]+dplanck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (1./absList[:-1]+1./absList[1:])/2.
rosselandabs_list.append(dplanck_sum / numpy.sum(xx))

r = 20.
xList = 2.*Pi *r/waterice[:,0]
for i in range(len(xList)):
    refrel = complex(refRe[i],refIm[i])
    qabs,qsca,g = bhmie(xList[i],refrel)
    absList[i] = qabs
    scaList[i] = qsca
    gList[i]   = g
ax.semilogx(waterice[:,0], absList, 'k-', label=r'r=20$\mu$m')
xx = (planck_list[:-1]+planck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (absList[:-1]+absList[1:])/2.
absmean_list.append(numpy.sum(xx) /planck_sum)
xx = (dplanck_list[:-1]+dplanck_list[1:])/2.*(wn[:-1] - wn[1:]) * \
    (1./absList[:-1]+1./absList[1:])/2.
rosselandabs_list.append(dplanck_sum / numpy.sum(xx))

print absmean_list, rosselandabs_list

ax.set_ylim([1e-5, 10])
ax.set_yscale('log')
ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel('Absorption efficiency')
ax.grid()
ax.tick_params(which='both',direction='out')
ax.legend(loc=0,fontsize=10)
# ax = fig.add_subplot(312)
# ax.semilogx(waterice[:,0], gList)
# ax.set_xlabel('Wavelength (micron)')
# ax.set_ylabel('Asymmetry factor')

# ax = fig.add_subplot(313)
# ax.semilogx(waterice[:,0], scaList/extList)
# ax.set_xlabel('Wavelength (micron)')
# ax.set_ylabel('Single scattering albedo')

# fig, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)
# # fig.set_size_inches(12, 6)
# cmap = mpl.cm.spectral #RdBu_r 
# levelp = np.logspace(2, np.log10(1e5), 40)

# CS2 = ax.contourf(X,Y, np.mean(y1out,axis=2) /86400,  
#     # reff*1e6, levels=np.linspace(5,40, 36), 
#     locator=mpl.ticker.LogLocator(), levels=np.logspace(-1,3,4*4+1), 
#     cmap=cmap  )
# ax.contour(X,Y, np.mean(y1out,axis=2) /86400, levels=[tau_auto[i] / 86400], colors='b',
#     linewidths=2, linestyles='--')
# ax.set_yscale('log')
# ax.set_ylim([1e5, 4e2])
# ax.set_xlim([-90, 90])
# if (i==0):
#     ax.set_ylabel(r'Pressure (Pa)')
# ax.tick_params(which='both',direction='out')
# ax.grid()
# # CB = fig.colorbar(CS2, ax=ax, orientation='horizontal')
# # ax.text(-90, 3e2, titles[i], fontsize=12)
# f.close()
# CB = fig.colorbar(CS2, ax=axes[:].ravel().tolist(),location='right')#orientation='horizontal')
# CB.ax.set_ylabel(r'Time$_{grav}$ (day)')
# mass_to_depth = 3./4*1e2 #from cloud mass to optical thickness
# CB = fig.colorbar(CS2, ax=axes.ravel().tolist(), pad=0.1, location='bottom')
# CB.set_ticks()
# print (3*2e-15 /(4*np.pi*CCN *1e3))**(1./3) *1e6
# plt.tight_layout()
# fig.savefig('CLD_MASS_EQ.pdf')

plt.show()
module manabe_convection_mod

!19/08/05, dilute form

use simple_sat_vapor_pres_mod, only: escomp, descomp

use constants_mod, only: grav, rdgas, rvgas, kappa, latent_heat, &
        cp_air, cp_vapor, cp_water, &
        HLV

implicit none
private

public :: qsat, compute_k, warm_start, moist_convection

!real(8) :: tau_bm = 7200.

real(8), public, parameter :: d622 = rdgas/rvgas
real(8), public, parameter :: d378 = 1.-d622
real(8), public, parameter :: d608 = d378/d622
!##==============================================
!##==============================================
contains

subroutine compute_geopotential(t, q, p_half, p_full, z_full, z_half)

real(8), intent(in), dimension(:) :: t, q, p_full
real(8), intent(in), dimension(:) :: p_half
real(8), intent(out), dimension(:) :: z_full
real(8), intent(out), dimension(:) :: z_half

integer :: n, k
real(8) :: virtual_t

n = size(t)
z_half(n+1) = 0.
    !virtual_t = t * (1. + d608 * q)
    !for k in range(len(z_full) - 1, 0, -1):
    do k = n, 2, -1
        virtual_t = t(k) * (1. + d608 * q(k)) 
        z_half(k) = z_half(k+1) + rdgas* virtual_t *&
                    (log(p_half(k+1) / p_half(k)))

        z_full(k) = z_half(k+1) + rdgas* virtual_t *&
                    (log(p_half(k+1) / p_full(k)))
    end do

    virtual_t = t(1) * (1. + d608 * q(1)) 
    z_full(1) = z_half(2) + rdgas* virtual_t *&
                (log(p_half(2) / p_full(1)))
    z_half(1) = 0.
    !return z_full, z_half
end subroutine compute_geopotential
!##==============================================
!hc = 1.
!do_evap = False

subroutine qsat(T,p, qs)
 
real(8), intent(in) :: T, p
real(8), intent(out) :: qs
real(8) :: esat
    call escomp(T, esat)
    qs = 1.
    if (esat .le. p) qs = d622 * esat / (p - d378 * esat)
    !return qs
end subroutine qsat

subroutine dqsatdT(T,p,dqs)

real(8), intent(in) :: T, p
real(8), intent(out):: dqs
real(8) :: esat, desat
    call escomp(T, esat)
    call descomp(T, desat)
    dqs = 1.e6
    if ( esat .le. p) dqs = d622 * p / (p - d378 * esat)**2 * desat
end subroutine dqsatdT

!#considering liquid in the air
subroutine rsat(T,p, rs)

real(8), intent(in) :: T, p
real(8), intent(out):: rs
real(8) :: esat
    call escomp(T, esat)
    rs = 0. !-1e5
    if (esat .le. p) rs = d622 * esat / (p -  esat)
end subroutine rsat

subroutine drsatdT(T,p, drs)

real(8), intent(in) :: T, p
real(8), intent(out):: drs
real(8) :: esat, desat
    call escomp(T, esat)
    call descomp(T, desat)
    if ( esat .le. p) drs = d622 * p / (p - esat)**2 * desat
end subroutine drsatdT

subroutine ff(dT, dp1, dp2, p1, p2, Dry_m, kini, x, val, lmass1,lmass2)
real(8), intent(in) :: dT, dp1, dp2, p1, p2, Dry_m, kini, x
real(8), intent(out):: val, lmass1, lmass2
real(8) :: t1n, t2n, yy, r1n, r2n, lh, rt1, rt2, &
       dmass1, dmass2, kfin 

t2n = x - dT
t1n = x !t2n + dT
!call rsat(t1n, p1, r1n)
!call rsat(t2n, p2, r2n)
!yy = (dp1/(1.+r1n) + dp2/(1.+r2n))/Dry_m -1.
!rt1 = r1n + yy * (1.+r1n)
!rt2 = r2n + yy * (1.+r2n)
!dmass1 = dp1 / (1.+rt1) !*(rt1 - r1n) 
!dmass2 = dp2 / (1.+rt2) !*(rt2 - r2n) 
!lmass1 = dp1 / (1.+rt1) *(rt1 - r1n) 
!lmass2 = dp2 / (1.+rt2) *(rt2 - r2n) 
!final moist enthalpy
!kfin = cp_air *  (dmass2*t2n+dmass1*t1n) + &
!        cp_water*(dmass2*rt2*t2n + dmass1*rt1*t1n)
!call latent_heat(t2n, lh)
!kfin = kfin + lh*dmass2*r2n
!call latent_heat(t1n, lh)
!kfin = kfin + lh*dmass1*r1n
!print *, kfin

call qsat(t1n, p1, r1n)
call qsat(t2n, p2, r2n)
kfin = cp_air*(t1n*dp1+t2n*dp2) + hlv*(r1n*dp1+r2n*dp2)
lmass1=Dry_m - (r1n*dp1+r2n*dp2)
lmass2=0.

val = kfin - kini 

end subroutine ff

subroutine newtSolve_ff(dT, dp1, dp2, &
                      q1, q2, T1, T2, p1, p2, T1new, T2new, lmass1, lmass2)
real(8), intent(in) :: dT, dp1, dp2, q1, q2, T1, T2, &
        p1, p2
real(8), intent(out):: T1new, T2new, lmass1, lmass2
real(8) :: x, val, x1, val1, x2, val2, deriv
real(8) :: Dry_m, kini, lh
real(8) :: epss = 1.e-7
integer :: nmax = 100, n 

!Dry air mass in the two layers
!Dry_m = dp1*(1.-q1) + dp2*(1. - q2)
!initial moist enthalpy
!kini = cp_air *  (dp2*(1-q2)*T2+dp1*(1-q1)*T1) + &
!        cp_water*(dp2*q2*T2 + dp1*q1*T1)
!call latent_heat(T2, lh)
!kini = kini + lh*dp2*q2
!call latent_heat(T1, lh)
!kini = kini + lh*dp1*q1
Dry_m = dp1*q1 + dp2*q2
kini = cp_air*(T1*dp1+T2*dp2) + hlv*(q1*dp1+q2*dp2)

    x = T1 !2
    x1 = x - 0.1
    val = 10.
    n = 0
    do while (abs(val) .gt. 1.e-3)
       call ff(dT, dp1, dp2, p1, p2, Dry_m, kini, x, val, lmass1,lmass2)
       
       x1 = x + epss
       call ff(dT, dp1, dp2, p1, p2, Dry_m, kini, x1, val1, lmass1,lmass2)
       x2 = x - epss
       call ff(dT, dp1, dp2, p1, p2, Dry_m, kini, x2, val2, lmass1,lmass2)
       deriv = (val1 - val2) / 2./ epss
       !deriv = (val1 - val) / (x1 - x)
       !x1 = x
       x = x - val / deriv  
       !print *, val, deriv, x
       
       n = n+1
       if (n > nmax) then
          print *, 'does not converge', x, lmass1
          lmass1 = -1e5
          exit
       end if  
    end do
    
    T2new = x - dT
    T1new = x !+ dT
    call ff(dT, dp1, dp2, p1, p2, Dry_m, kini, x, val, lmass1, lmass2)
    !print *, 'final: ', val, deriv, x

    !if (isnan(lmass1)) then
    !   print *, 'rain NaN'
    !   lmass1 = -1e5
    !endif

    call rsat(T1new, p1, val1)
    if (val1 .le. 0) lmass1 = -1e5
end subroutine newtSolve_ff

!=====================================================================

subroutine compute_k(tin,qin,phalf, k)

real(8), intent(in), dimension(:) :: tin, qin
real(8), intent(in), dimension(:) :: phalf
real(8), intent(out) :: k
real(8) :: kpum
integer :: i
    k = 0.
    do i = 1, size(tin)
       kpum = (cp_air*(1 - qin(i)) + cp_vapor*qin(i))*tin(i)
       k = k + kpum * (phalf(i+1) - phalf(i)) / grav  
    end do
    !return sum(k* (phalf[1:] - phalf[:-1])/cons.grav)
end subroutine compute_k

!#initial profile
subroutine dlnpdlnT_initial(lnT,lnp, slope)

real(8), intent(in) :: lnT, lnp
real(8), intent(out):: slope
real(8) :: T, p, qs, D, rs, pa, es, A, B, C, ll

    T = exp(lnT)
    p = exp(lnp)
    
    call qsat(T,p, qs)
    !call latent_heat(T, ll)
    ll = hlv *1.
    D = ll / (rvgas * T)

    if (qs .eq. 1) then
        slope = D
    else
        rs = qs / (1. - qs)
        pa = p / (rs /d622 + 1.)
        es = p - pa
        A = ll * rs / rdgas / T
        B = rs * (cp_vapor / cp_air + &
              (ll / rvgas / T-1.) * ll / cp_air / T)
        C = 1. / (kappa * (1. + A) / (1. + B)) !#dlnpa / dlnT
        
        slope = (pa * C + es * D) / p
    end if
end subroutine dlnpdlnT_initial


subroutine rk4_f(x0, y0, dx, nx, x, y)

real(8), intent(in) :: x0, y0, dx
integer, intent(in) :: nx
real(8), intent(out), dimension(nx+1) :: x, y
integer :: n
real(8) :: k1, k2, k3, k4

x(1) = x0
y(1) = y0

do n = 1, nx
   call dlnpdlnT_initial(x(n), y(n), k1)
   call dlnpdlnT_initial(x(n) + dx/2., y(n) + dx/2.*k1, k2)
   call dlnpdlnT_initial(x(n) + dx/2., y(n) + dx/2.*k2, k3)
   call dlnpdlnT_initial(x(n) + dx, y(n) + dx *k3, k4)
   x(n+1) = x(n) + dx
   y(n+1) = y(n) + dx/6.*(k1+k4 + 2.*(k2+k3))
end do
end subroutine rk4_f

subroutine warm_start(ts, pd, nx, psg, tt, lnpp, qs)
real(8), intent(in) :: ts, pd
integer, intent(in) :: nx
real(8), intent(out):: psg
real(8), intent(out), dimension(nx+1) :: tt, lnpp, qs
integer :: i
real(8) :: psc, x1, x2, dx
real(8), dimension(nx+1) :: lntlist, lnplist, lntt, pp

    call escomp(ts, psc)
    psg = pd + psc !sat.escomp(ts)

    x1 = log(ts)
    x2 = log(20.)
    !nx = 400
    dx = (x2 - x1) /nx
    call rk4_f(x1, log(psg), dx, nx, lntlist, lnplist)
    do i = 1, nx+1
       lntt(i) = lntlist(nx+2-i)
       lnpp(i) = lnplist(nx+2-i)
    end do
    tt = exp(lntt)
    pp = exp(lnpp)
    do i = 1, nx+1
       call qsat(tt(i), pp(i), qs(i))
    end do
    !lntlist = [x1]
    !lnplist = [np.log(psg)]
    !int_g = integrator(dlnpdlnT_initial, x1, np.log(psg), dx)
    !for i in range(nx):
        !lnt, lnp = int_g.next()
        !lntlist.append(lnt)
        !lnplist.append(lnp)
    !lntlist.reverse()
    !lnplist.reverse()
    !lntt = np.array(lntlist)
    !lnpp = np.array(lnplist)
    !tt = np.exp(lntt)
    !pp = np.exp(lnpp)
    !qs = tt *1.
    !for i in range(len(tt)):
    !    qs[i] = cond.qsat(tt[i], pp[i])

    !return psg, tt, lnpp, qs
end subroutine warm_start

!#==================================================================    
subroutine dlnpdlnT(T,p, slope)

real(8), intent(in) :: T, p
real(8), intent(out):: slope
real(8) :: qs, lh, D, rs, pa, es, A, B, C

    call qsat(T,p, qs)
    !call latent_heat(T, lh)
    lh = hlv *1.
    D = lh / (rvgas * T)

    if (qs .eq. 1) then
       slope = D
    else
        rs = qs / (1. - qs)
        pa = p / (rs /d622 + 1.)
        es = p - pa
        A = lh * rs / rdgas / T
        B = rs * (cp_vapor / cp_air + &
              (lh / rvgas / T-1.) * lh / cp_air / T)
        C = 1. / (kappa * (1. + A) / (1. + B)) !#dlnpa / dlnT
        slope =  (pa * C + es * D) / p
    end if
end subroutine dlnpdlnT

subroutine convec(T1,P1,P2,PH, T2)

real(8), intent(in) :: T1, P1, P2, PH
real(8), intent(out):: T2
real(8) :: lr1, TH, lr2
    call dlnpdlnT(T1,P1, lr1)

    TH = T1 * (PH/P1)**(1./lr1)

    call dlnpdlnT(TH,PH, lr2)

    T2 = T1 * (P2 / P1)**(1./lr2)
end subroutine convec
    !return T2
!##==============================================
subroutine moist_convection(tin, qin, p_full, p_half, &
                rain, Tref, qref, conv_cloud) !, p_full, p_half, Ep, rain_profile)

real(8), intent(in), dimension(:,:,:) :: tin, qin, p_full
real(8), intent(in), dimension(:,:,:) :: p_half
real(8), intent(out),dimension(:,:)   :: rain!, Ep
real(8), intent(out),dimension(:,:,:) :: Tref, qref, conv_cloud!, p_full, rain_profile
!real(8), intent(out), dimension(:) :: p_half
!real(8) :: k1, kk1, k2, k3
real(8) :: T1, P1, P2, PH, T2, T2d
real(8) :: q1, q2, deltap1, deltap2, fq
real(8) :: cp1,cp2, cp3, cp4, deltat
real(8) ::  DTref, T1new, T2new
real(8) :: dp, lmass, tcon, dpn, kappa2
real(8) :: lh, lmass1, lmass2
!real(8), allocatable, dimension(:) :: pfulln, kappa1, zf
!real(8), allocatable, dimension(:) :: phalfn, zh
integer :: iter1, k, i, n, j, nc

    n = size(tin)
    !allocate(pfulln(n))
    !allocate(kappa1(n))
    !allocate(zf(n))
    !allocate(phalfn(n+1))
    !allocate(zh(n+1))
    !Ep = 0.
    rain = 0.
!    Ep_2 = 0. !check energy conservation
    Tref = tin *1.
    qref = qin *1.
    conv_cloud = 0.
    !p_full = p_full_in *1.
    !p_half = p_half_in *1.
    !rain_profile = 0.
        
    do iter1 = 1, 1 !20 !5
       do i = 1, size(p_full,1)
       do j = 1, size(p_full,2)
       do k = size(p_full,3), 2, -1
          T1 = Tref(i,j,k)
          P1 = p_full(i,j,k)
          P2 = p_full(i,j,k-1)
          PH = p_half(i,j,k)
          deltap1 = p_half(i,j,k+1) - p_half(i,j,k)
          deltap2 = p_half(i,j,k) - p_half(i,j,k-1)

          T2d = T1 * (P2 / P1)**kappa
          dry_conv: if (Tref(i,j,k-1) .lt. T2d) then
                   DTref = T1 - T2d
                   deltat = (Tref(i,j,k-1) - Tref(i,j,k) + DTref)*deltap2 / &
                                (deltap1 + deltap2)
                   Tref(i,j,k) = T1 + deltat
                   Tref(i,j,k-1) = T2d + deltat
                   fq = (qref(i,j,k-1)*deltap2 + qref(i,j,k)*deltap1) / &
                                (deltap1 + deltap2)
                   !qref(i,j,k)   = fq
                   !qref(i,j,k-1) = fq 

          else
          call convec(T1, P1, P2, PH, T2)
          !call compute_k(Tref, qref, p_half, k1)

          lapserate: if (Tref(i,j,k-1) .lt. T2) then
             DTref = T1 - T2
             
             call newtSolve_ff(DTref, deltap1, deltap2, &
                      qref(i,j,k), qref(i,j,k-1), Tref(i,j,k), Tref(i,j,k-1), &
                      p_full(i,j,k), p_full(i,j,k-1),T1new, T2new,lmass1,lmass2)

             !negative rain->do nothing
             pos_rain: if (lmass1 .gt. 0) then 

                if (T1new .gt. 400.) print *, "too warm"
                Tref(i,j,k) = T1new
                Tref(i,j,k-1) = T2new
                call qsat(Tref(i,j,k), p_full(i,j,k), qref(i,j,k))
                call qsat(Tref(i,j,k-1), p_full(i,j,k-1), qref(i,j,k-1)) !new qsat after removal
                rain(i,j) = rain(i,j) + lmass1 / grav

                conv_cloud(i,j,k) = conv_cloud(i,j,k) + lmass1 / (deltap1+deltap2)
                conv_cloud(i,j,k-1)=conv_cloud(i,j,k-1)+ lmass1/(deltap1+deltap2)
          
             end if pos_rain!moist convection
          end if lapserate

          end if dry_conv
          !if (isnan(Tref(i,j,k))) print *, k, rain(i,j)     
  
       end do
       end do
       end do                
    end do

end subroutine moist_convection

 
subroutine cold_trap(t, q, phalf, pfull, nc)
real(8), intent(inout), dimension(:) :: t, q
real(8), intent(in), dimension(:) :: phalf
real(8), intent(in), dimension(:) :: pfull
integer, intent(in) :: nc

real(8) :: q2
integer :: i

do i = nc-1, 1, -1
   call qsat(t(i), pfull(i), q2)
   if (q2 .gt. q(i+1)) exit
   q(i) = q2
end do

q(1:i) = q(i+1)
end subroutine cold_trap

                    
    !return rain, Tref, qref, p_full, p_half, Ep, Ep_2 , rain_profile

end module manabe_convection_mod


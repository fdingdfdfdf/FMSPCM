module cloud_evap_mod

!20/10/13, DF           

use simple_sat_vapor_pres_mod, only: escomp, descomp

use constants_mod, only: grav, rdgas, rvgas, kappa, latent_heat, &
        cp_air, cp_vapor, cp_water, &
        HLS, HLV, TFREEZE, r_cld, dens_h2o

implicit none
private

public :: calculate_rehum, cloud_evap

real(8), public, parameter :: d622 = rdgas/rvgas
real(8), public, parameter :: d378 = 1.-d622
real(8), public, parameter :: d608 = d378/d622
real(8), parameter :: ka0 = 2.43e-2 !J/(m s K) thermal conductivity
real(8), parameter :: psi0 = 2.26e-5 !m**2/s diffusivity
real(8), parameter :: mu0 = 1.8325 !viscosity
real(8), parameter :: mu1 = 1.718
real(8), parameter :: tt1 = 296.16
real(8), parameter :: cc  = 120
!##==============================================
!##==============================================
contains

subroutine calculate_rehum(t, q, p_full, rehum)
real(8), intent(in), dimension(:,:,:) :: t, q, p_full
real(8), intent(out), dimension(:,:,:) :: rehum

real(8), dimension(size(t,1),size(t,2),size(t,3)) :: esat, ewv

call escomp(t, esat)
ewv   = q*p_full/(d622+d378*q)
rehum = ewv / esat

end subroutine calculate_rehum
!##==============================================
subroutine cloud_evap(dt, q_cld, t, q, p_full, &
        rehum, cld_dt_tg, cld_dt_qg)
real(8), intent(in) :: dt
real(8), intent(in), dimension(:,:,:)  :: q_cld, t, q, p_full, rehum
real(8), intent(out), dimension(:,:,:) :: cld_dt_tg, cld_dt_qg
real(8), dimension(size(t,1),size(t,2),size(t,3)) :: esat, denom

! eq. 37, wisner et al. 1972
call escomp(t, esat)
!where (t .gt. TFREEZE)
!   lh = HLV
!elsewhere
!   lh = HLS
!endwhere
denom = HLV**2/rvgas/t**2/ &
        (ka0*mu0/mu1*(tt1+cc)/(t+cc)*(t/tt1)**1.5) + &
        rvgas*t/esat/(psi0*1e5/p_full*(t/TFREEZE)**1.81)
cld_dt_qg = q_cld *3 *(1-rehum) /dens_h2o /r_cld**2 /denom *dt

where(cld_dt_qg .gt. q_cld) cld_dt_qg = q_cld
where(cld_dt_qg .lt. 0)     cld_dt_qg = 0.
cld_dt_tg = 0.-HLV* cld_dt_qg /cp_air !evaporation cools the air

end subroutine cloud_evap

end module cloud_evap_mod


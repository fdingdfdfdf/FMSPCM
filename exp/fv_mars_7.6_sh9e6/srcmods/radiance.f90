module radiance_mod
! ==================================================================================
! index order is different from the main program
! 1-- surface now
! nLay-- TOA
! ==================================================================================

   use            fms_mod, only: set_domain, write_version_number, &
                              mpp_pe, mpp_root_pe, error_mesg, FATAL, WARNING

   use            fms_mod, only: stdlog, check_nml_error, close_file,&
                              open_namelist_file, stdout, file_exist, &
                              read_data, write_data, open_file, &
                              nullify_domain, lowercase

   use constants_mod,      only: stefan, cp_air, grav, pstd_mks, &
                                B_nu, rho_cp, &
                                c, c2, atm, Rstar, Tref, mpr, &
                                nS, datadir, kappa, &
                                dens_h2o, tfreeze, r_cld, &
                                seconds_per_day, deltat_rad, sigma_sca

   use line_profs,         only: lorentz, doppler, voigt
   use cia,                only: setup_cia, calculate_CIA, &
                                calculate_CO2CO2_CIA, &
                                interpolateH2Ocont_PPC, interpolateH2Ocont_CKD, &
                                calculate_O2O2_cia
   use mult_scat_core, only : get_Iup_Idn, setup_mult_scat_core
                              !backscatter, backscatter_table , &
                              !angular_backscatter 

   use fv_mp_mod,             only: domain, gid, masterproc, &
                             mp_reduce_sum, mp_reduce_min, mp_reduce_max
   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)
!==================================================================================
implicit none
private
!==================================================================================
! version information 
character(len=128) :: version='$Id: radiance.f90 $'
character(len=128) :: tag='homemade'
!==================================================================================
! public interfaces
public :: radiance_init , radiance_update  ! 
!==================================================================================
  real    :: Trot   = 1. !Earth day
  real    :: Torb   = 365. !Earth day
  real    :: starlum= 3.828e26
  real    :: semimajor = 1.496e11 !meter
  real    :: ecc    = 0.
  real    :: obl    = 0. !deg
  real    :: prec   = 180.  !deg
  real    :: dist, sublat, sublon, kappa_ang, kappa_new, sublon_new

  real    :: A_ice   = 0.5
  real    :: A_water = 0.05 !ocean albedo
  logical :: adleo = .false.
  real    :: gsca !asymmetry factor
!  logical :: tidal_lock= .true.
  integer :: nGas      = 3           ! number of species

  integer :: nAng      = 1 !4            ! number of angles to integrate over
  integer :: nlinesMAX = 3000000      ! maximum number of lines to read
  integer :: nTem      = 1 !5            ! number of cross-section temperatures at each layer (1 or 2)
  
  real    :: solar_constant  = 1093. !1360.0

  logical :: master
! module variables
  logical :: initialized =.false.
  integer :: vgas = 1 !0, dry 
!  real    :: noon_longitude  = 270.
!  real    :: del_sol         = 1.4
! modif omp: winter/summer hemisphere
!  real    :: del_sw          = 0.
  real    :: kappa_gray = 1e-4 !.5d-4 !0.01 !5.0d-3 

  integer  :: nlines_beguier = 7786 ! number of lines in Beguier+ JQRST (2015) CH4 data
  logical  :: use_beguier_2015 = .true. ! include Beguier+ JQRST (2015) CH4 data in near-IR?
  logical  :: HITEMP = .False.                  ! use HITEMP line lists?
!  real :: deltanu_trunc = 25.0d0      ! truncation wavenumber [cm^-1]
  real     :: dTlay = 25.0 !50.0 !25.0
  real, parameter :: Omega_stel = 6.8d-5  ! stellar solid angle [sr]
  ! Sun at Earth by default, e.g.: http://www.geo.mtu.edu/~scarn/teaching/GE4250/radiation_lecture_slides.pdf
  logical :: stellar_blackbody = .false.  ! use a blackbody curve for the stellar spectrum?
  integer, parameter :: nS_stel = 50000             ! number of points in stellar spectrum

!  integer :: iGas_H2O = 1
!  integer :: iGas_CO2 = 2
!  integer :: iGas_O2  = 3
!  integer :: iGas_CH4 = 6
!  integer :: iGas_N2  = 7
!  integer :: iGas_O3  = 8
    !------- inputs --------
  real :: nu_lw1 = 1.0 !1.0
  real :: nu_lw2 = 3000.0 !5000.0
  logical :: gray_debug   = .false. ! uniform opacity vs. wavenumber
  !------- namelist parameters --------
  real :: nu_sw1 = 1e3 !61.0 !1.0                               ! longwave starting wavenumber [cm^-1]
  real :: nu_sw2 = 5e4 !30000.0 !50000.0                              ! longwave finishing wavenumber [cm^-1]
  real :: Asurf  = 0.20                             ! surface albedo []
!  real :: Fstel0                               ! stellar flux [W/m2]
!  real :: cosa0                                ! stellar zenith angle cosine []
  logical ::  rayleigh_top                         ! apply Rayleigh scattering at TOA (increases albedo)?
  logical ::  multiple_scat
    ! constants once loaded
  integer, allocatable, dimension(:) :: nlines, deltanu_trunc            ! actual number of lines (nGas
  real,    allocatable, dimension(:)    :: mu_i !nGas
  integer, allocatable, dimension(:,:) ::  mol,  iso    ! molecule number (nGas, nlinesMAX)
  real,    allocatable, dimension(:,:) ::  nu0, Sref, einstein, gam_air, gam_self,&
          Egnd, ncoeff, delta , &
          gam_lore, gam_dopp, Strue !(nGas, nlinesMAX)

!  integer iLay, iLev, iAng, iS                    ! do loop variables
!  integer iRev, iGas, ierr

  real, allocatable, dimension(:) :: cosa, ang_wt !nAng
  real, allocatable, dimension(:) ::  nu_lw, nu_sw, sigma_lw, & !nS  
          sigma_CIA, sigma_sw
  real :: dnu_lw, dnu_sw

  !vertical profile for rad
  real, allocatable, dimension(:,:,:) ::  play, Tlay, lnplay
  real, allocatable, dimension(:,:,:) ::  Tlev, plev, lnplev 
  real, allocatable, dimension(:,:) :: t_surf, dt_surf_rad, &
       net_surf_sw_down, surf_lw_down 
  real, allocatable, dimension(:,:) :: Fstel0, cosa0
!  real, allocatable, dimension(:,:) :: ss, cos_lat, p2, solar
  real, allocatable, dimension(:,:) :: T_k_grid !nlev, nTem
  integer, allocatable, dimension(:,:,:) :: iT1, iT2, iT1_out
  real, allocatable, dimension(:,:,:,:) ::  T_lin_weight !(nS,nLay) 
  real, allocatable, dimension(:,:,:,:) :: f_i !nGas
  real, allocatable, dimension(:,:,:)   :: dp, mu_avg
  real, allocatable, dimension(:,:,:,:) :: sigma_lw_ar, sigma_sw_ar 
  real, allocatable, dimension(:,:,:)   :: sigma_total_dry, dtau_sw, & 
        sigma_sw_total_dry
  real, allocatable, dimension(:,:,:,:) :: sigma_t_i, &
          dtau_lw, dtau_lw_a, dTran, &
          Bnu_lay, Cterm, sigma_d_i, sigma_v_i, &
          dtau_sw_i, dTran0, &
          dtau_sw_cld, dtau_sw_ray
  real, allocatable, dimension(:,:,:,:) :: tau_lw, Bnu_lev, &
          I_lev_up, I_lev_dn, I_lev_dir
  real, allocatable, dimension(:,:,:,:,:) :: log_sig_d, log_sig_v,log_dtau_sw 

real, allocatable, dimension(:,:,:) :: OLRnu, OLRnu_temp, OSRnu, OSRnu_temp, GSRnu, GSRnu_temp, Bnu_s, tau_lw_inf, &
Beff, I_stel_dn, F_stel_dn, tau_sw_inf, Rs, Rs_out, &
Rpla, Rpla_out
real, allocatable, dimension(:,:) :: OLR_temp, OSR_temp, GSR_temp, OLR, OSR, GSR
real, allocatable, dimension(:,:,:) :: Ilev_up_lw, Ilev_dn_lw, Flev_up_lw, Flev_dn_lw, &
       Ilev_up_sw, Ilev_dn_sw, Flev_up_sw, Flev_dn_sw , Ilev_dir_sw, Flev_dir_sw
real, allocatable, dimension(:,:,:) :: dF_lw, dF_sw, dTdt_lw, dTdt_sw, dTdt

    ! parameters for multiple scattering calculation 
    real,allocatable,dimension(:,:,:) ::  aL, bL, b0L   ! single scattering albedo [],backscatter parameter [] after delta-two-term rescaling
    real,allocatable,dimension(:,:,:) :: ssa_sw_i, delta_f, x1_rescaled, f_cloud
    !real,allocatable,dimension(:,:) :: b0L_cld
    real :: bL_cld       ! asymmetry parameter []
    real,allocatable,dimension(:,:,:) ::  tauL,  I_lev_up_ms, I_lev_dn_ms, I_lev_dir_ms ! optical depth (0 at TOA, max. at BOA) []

real, save :: pi, deg_to_rad , rad_to_deg

real :: mu_dry
   real, allocatable, dimension(:) :: vgas_mask
   character(len=3),  allocatable, dimension(:) :: gas_name
   real(8),  allocatable, dimension(:)  :: gas_molarconc
   character(len=3), dimension(10)      :: gas_name_MAX
   real(8),          dimension(10)      :: gas_molarconc_MAX,deltanu_trunc_MAX
   ! initialize all these variables to -1
   integer :: iGas_H2   = -1
   integer :: iGas_He   = -1
   integer :: iGas_H2O  = -1
   integer :: iGas_CO2  = -1
   integer :: iGas_CO   = -1
   integer :: iGas_N2   = -1
   integer :: iGas_O2   = -1
   integer :: iGas_O3   = -1
   integer :: iGas_SO2  = -1
   integer :: iGas_H2S  = -1
   integer :: iGas_CH4  = -1
   integer :: iGas_NH3  = -1
 
 namelist/radiance_nml/  nGas, nAng, nTem, &
         gas_name_MAX, gas_molarconc_MAX, deltanu_trunc_MAX, &
         stellar_blackbody, gray_debug, rayleigh_top, multiple_scat, &
         Asurf, A_ice, A_water, gsca, adleo, Trot, Torb, semimajor, &
         starlum, ecc, obl, prec
!==================================================================================
!-------------------- diagnostics fields -------------------------------
integer :: id_flux_lw, id_tau_lw_inf, id_flux_sw, &
        id_sigma_lw, id_sigma_di, id_sigma_vi
integer :: id_OLRnu, id_OLR !, id_solar, id_OSR, id_OSRnu 
integer :: id_OSRnu, id_OSR, id_solar, id_solarnu, id_cosa
integer :: id_hrsw, id_hrlw
integer :: id_Fdn_sw, id_Fup_sw, id_Fdir_sw
logical :: used
character(len=10), parameter :: mod_name = 'lbl'

real :: missing_value = -999.
contains
! ==================================================================================
! ==================================================================================

subroutine radiance_init(is, ie, js, je, num_levels,&
        axes, Time, agrid, q, pt, ps, pfull)
!-------------------------------------------------------------------------------------
integer, intent(in), dimension(5) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real, intent(in), dimension(:,:,:) ::   pt, q, pfull
real, intent(in), dimension(:,:,:) ::   agrid
real, intent(in), dimension(:,:)   ::   ps
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
logical :: water_file_exists
integer :: i,j, k, iGas, countgas
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile
unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=radiance_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'radiance_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=radiance_nml)
!unit = open_file ('input.nml', action='read')
!ierr=1
!do while (ierr /= 0)
!   read  (unit, nml=radiance_nml, iostat=io, end=10)
!   ierr = check_nml_error (io, 'radiance_nml')
!enddo
!10 call close_file (unit)

!unit = open_file ('logfile.out', action='append')
!if ( mpp_pe() == 0 ) then
!  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
!  write (unit, nml=radiance_nml)
!endif
!call close_file (unit)

    if((multiple_scat.eqv..true.) .and. (nAng.ne.1))then
       call error_mesg('radiance_mod:','Currently the code only runs with two-stream multiple',FATAL)    
       !print*,'scattering. Set nAng = 1 in params_consts.f90.'
       !stop
    endif

pi    = 4.0*atan(1.)
deg_to_rad = 2.*pi/360.
rad_to_deg = 360.0/2./pi

initialized = .true.
master = (gid == masterproc)
allocate(vgas_mask      (nGas))
allocate(nlines         (nGas))
allocate(deltanu_trunc  (nGas)) ! actual number of lines (nGas
allocate(mu_i           (nGas))
allocate(mol            (nGas, nlinesMAX))
allocate(iso            (nGas, nlinesMAX)) ! molecule number (nGas, nlinesMAX)
allocate(nu0            (nGas, nlinesMAX))
allocate(Sref           (nGas, nlinesMAX))
allocate(einstein       (nGas, nlinesMAX))
allocate(gam_air        (nGas, nlinesMAX))
allocate(gam_self       (nGas, nlinesMAX))
allocate(Egnd           (nGas, nlinesMAX))
allocate(ncoeff         (nGas, nlinesMAX))
allocate(delta          (nGas, nlinesMAX))
allocate(gam_lore       (nGas, nlinesMAX))
allocate(gam_dopp       (nGas, nlinesMAX))
allocate(Strue          (nGas, nlinesMAX))

allocate(cosa           (nAng))
allocate(ang_wt         (nAng))
allocate(nu_lw          (nS))
allocate(nu_sw          (nS))
allocate(sigma_lw       (nS))
allocate(sigma_sw       (nS))
allocate(sigma_CIA      (nS))

allocate(play           (ie-is+1, je-js+1, num_levels))
allocate(Tlay           (ie-is+1, je-js+1, num_levels))
allocate(lnplay         (ie-is+1, je-js+1, num_levels))
allocate(plev           (ie-is+1, je-js+1, num_levels+1))
allocate(Tlev           (ie-is+1, je-js+1, num_levels+1))
allocate(lnplev         (ie-is+1, je-js+1, num_levels+1))

allocate(t_surf         (ie-is+1, je-js+1))
allocate(dt_surf_rad    (ie-is+1, je-js+1))
allocate(net_surf_sw_down (ie-is+1, je-js+1))
allocate(surf_lw_down     (ie-is+1, je-js+1))

allocate(T_k_grid       (num_levels, nTem))
allocate(iT1		(ie-is+1, je-js+1, num_levels))
allocate(iT2            (ie-is+1, je-js+1, num_levels))
allocate(T_lin_weight   (ie-is+1, je-js+1, num_levels, nS))

allocate (sigma_sw_ar           (num_levels, nS, nGas, nTem))
allocate (sigma_lw_ar           (num_levels, nS, nGas, nTem))
allocate (sigma_total_dry       (num_levels, nS, nTem))
!allocate (dtau_sw               (num_levels, nS, nTem))
allocate (sigma_sw_total_dry    (num_levels, nS, nTem))

allocate (log_sig_d     (ie-is+1, je-js+1, num_levels, nS, 2))
allocate (log_sig_v     (ie-is+1, je-js+1, num_levels, nS, 2))
!allocate (log_dtau_sw   (ie-is+1, je-js+1, num_levels, nS, 2))

allocate (dtau_sw_cld   (ie-is+1, je-js+1, num_levels, nS)) 
allocate (dtau_sw_ray   (ie-is+1, je-js+1, num_levels, nS)) 
allocate (dtau_sw_i     (ie-is+1, je-js+1, num_levels, nS)) !vertical
allocate (sigma_t_i     (ie-is+1, je-js+1, num_levels, nS))
allocate (sigma_d_i     (ie-is+1, je-js+1, num_levels, nS))
allocate (sigma_v_i     (ie-is+1, je-js+1, num_levels, nS))
allocate (dtau_lw       (ie-is+1, je-js+1, num_levels, nS))
allocate (dtau_lw_a     (ie-is+1, je-js+1, num_levels, nS)) !angle
allocate (dTran         (ie-is+1, je-js+1, num_levels, nS))
allocate (dTran0        (ie-is+1, je-js+1, num_levels, nS))
allocate (Bnu_lay       (ie-is+1, je-js+1, num_levels, nS))
allocate (Cterm         (ie-is+1, je-js+1, num_levels, nS))

allocate (tau_lw        (ie-is+1, je-js+1, num_levels+1, nS))
allocate (Bnu_lev       (ie-is+1, je-js+1, num_levels+1, nS))
allocate (I_lev_up      (ie-is+1, je-js+1, num_levels+1, nS))
allocate (I_lev_dn      (ie-is+1, je-js+1, num_levels+1, nS))
allocate (I_lev_dir     (ie-is+1, je-js+1, num_levels+1, nS))

allocate (OLRnu         (ie-is+1, je-js+1, nS))
allocate (OLRnu_temp    (ie-is+1, je-js+1, nS))
allocate (OSRnu         (ie-is+1, je-js+1, nS))
allocate (OSRnu_temp    (ie-is+1, je-js+1, nS))
allocate (GSRnu         (ie-is+1, je-js+1, nS))
allocate (GSRnu_temp    (ie-is+1, je-js+1, nS))
allocate (Bnu_s         (ie-is+1, je-js+1, nS))
allocate (Beff          (ie-is+1, je-js+1, nS))
allocate (tau_lw_inf    (ie-is+1, je-js+1, nS))
!allocate (tau_sw_inf    (ie-is+1, je-js+1, nS))
allocate (I_stel_dn     (ie-is+1, je-js+1, nS))
allocate (F_stel_dn     (ie-is+1, je-js+1, nS))
allocate (Rs            (ie-is+1, je-js+1, nS))
allocate (Rs_out        (ie-is+1, je-js+1, nS))
allocate (Rpla          (ie-is+1, je-js+1, nS))
allocate (Rpla_out      (ie-is+1, je-js+1, nS))

allocate (OLR         (ie-is+1, je-js+1))
allocate (OLR_temp    (ie-is+1, je-js+1))
allocate (OSR         (ie-is+1, je-js+1))
allocate (OSR_temp    (ie-is+1, je-js+1))
allocate (GSR         (ie-is+1, je-js+1))
allocate (GSR_temp    (ie-is+1, je-js+1))

allocate (Ilev_up_lw    (ie-is+1, je-js+1, num_levels+1))
allocate (Ilev_dn_lw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_up_lw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_dn_lw    (ie-is+1, je-js+1, num_levels+1))

allocate (Ilev_up_sw    (ie-is+1, je-js+1, num_levels+1))
allocate (Ilev_dn_sw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_up_sw    (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_dn_sw    (ie-is+1, je-js+1, num_levels+1))
!allocate (Ilev_dir_sw   (ie-is+1, je-js+1, num_levels+1))
allocate (Flev_dir_sw   (ie-is+1, je-js+1, num_levels+1))

allocate (dF_lw    (ie-is+1, je-js+1, num_levels))
allocate (dF_sw    (ie-is+1, je-js+1, num_levels))
allocate (dTdt_lw  (ie-is+1, je-js+1, num_levels))
allocate (dTdt_sw  (ie-is+1, je-js+1, num_levels))
allocate (dTdt     (ie-is+1, je-js+1, num_levels))
!allocate (Tcold    (ie-is+1, je-js+1, num_levels))
allocate (dp       (ie-is+1, je-js+1, num_levels))
allocate (mu_avg   (ie-is+1, je-js+1, num_levels))
allocate (f_i      (ie-is+1, je-js+1, num_levels, nGas))

allocate (aL       (ie-is+1, je-js+1, num_levels))
allocate (bL       (ie-is+1, je-js+1, num_levels))
allocate (b0L      (ie-is+1, je-js+1, num_levels))
allocate (tauL        (ie-is+1, je-js+1, num_levels+1))
allocate (f_cloud  (ie-is+1, je-js+1, num_levels))  !cloud mass concentration
allocate (ssa_sw_i (ie-is+1, je-js+1, num_levels))
allocate (delta_f  (ie-is+1, je-js+1, num_levels))
allocate (x1_rescaled (ie-is+1, je-js+1, num_levels))

allocate (I_lev_up_ms (ie-is+1, je-js+1, num_levels+1))
allocate (I_lev_dn_ms (ie-is+1, je-js+1, num_levels+1))
allocate (I_lev_dir_ms(ie-is+1, je-js+1, num_levels+1))

!allocate (ss               (ie-is+1, je-js+1))
!allocate (cos_lat          (ie-is+1, je-js+1))
!allocate (solar            (ie-is+1, je-js+1))
!allocate (p2               (ie-is+1, je-js+1))

allocate (Fstel0           (ie-is+1, je-js+1))
allocate (cosa0            (ie-is+1, je-js+1))
allocate(gas_name(nGas))
allocate(gas_molarconc(nGas))
!deltanu_trunc(:) = 25.0 !50.0 !25.0d0
!composition.f90
    ! create gas name and molar concentration arrays
    !if (.not.allocated(gas_name))      allocate(gas_name(nGas))
    !if (.not.allocated(gas_molarconc)) allocate(gas_molarconc(nGas))
    gas_name(1:nGas)      = gas_name_MAX(1:nGas)
    gas_molarconc(1:nGas) = gas_molarconc_MAX(1:nGas)
    vgas = 0
    ! identify the variable gas
    do iGas = 1, nGas
       if(gas_molarconc(iGas)<0.0d0 )then
          if (master) write(6,*) 'Variable gas is ',gas_name(iGas)
          vgas = iGas
          !call get_vgas_params
          exit
       end if
    end do
    if(vgas==0 .and. master) write(6,*) 'No variable gas found!'
    ! assign the 'iGas_X' labels
    do iGas=1,nGas
       if (gas_name(iGas)=="H2_") then
          iGas_H2=iGas
          mu_i(iGas_H2)=2.0159       ! molar mass of H2 [g/mol]
       elseif (gas_name(iGas)=="He_") then
          iGas_He=iGas
          mu_i(iGas_He)=4
       elseif (gas_name(iGas)=="H2O") then
          iGas_H2O=iGas
          mu_i(iGas_H2O)=18.01528     ! molar mass of H2O [g/mol]
       elseif (gas_name(iGas)=="CO2") then
          iGas_CO2=iGas
          mu_i(iGas_CO2)=44.0095      ! molar mass of CO2 [g/mol]
       elseif (gas_name(iGas)=="CO_") then
          iGas_CO=iGas
          mu_i(iGas_CO)=28
       elseif (gas_name(iGas)=="N2_") then
          iGas_N2=iGas
          mu_i(iGas_N2)=28.014       ! molar mass of N2 [g/mol]
       elseif (gas_name(iGas)=="O2_") then
          iGas_O2=iGas
          mu_i(iGas_O2)=32
       elseif (gas_name(iGas)=="O3_") then
          iGas_O3=iGas
          mu_i(iGas_O3)= 48.00        ! molar mass of O3 [g/mol]
       elseif (gas_name(iGas)=="SO2") then
          iGas_SO2=iGas
          mu_i(iGas_SO2)=64.06        ! molar mass of SO2 [g/mol]
       elseif (gas_name(iGas)=="H2S") then
          iGas_H2S=iGas
          mu_i(iGas_H2S)=34.08        ! molar mass of H2S [g/mol]
       elseif (gas_name(iGas)=="CH4") then
          iGas_CH4=iGas
          mu_i(iGas_CH4)=16.04        ! molar mass of CH4 [g/mol]
       elseif (gas_name(iGas)=="NH3") then
          iGas_NH3=iGas
          mu_i(iGas_NH3)=17
       endif
       f_i(:,:,:,iGas)=gas_molarconc(iGas)
       deltanu_trunc(iGas)=deltanu_trunc_MAX(iGas)
    enddo
!print *, gas_molarconc_MAX(iGas_CO2) !f_i(1,1,:,iGas_CO2)
!-----------------------------------------------------------------------
vgas_mask(1:nGas)=1.0
if(vgas>0) vgas_mask(vgas)=0.0
mu_dry = sum(gas_molarconc(:)*mu_i(:)*vgas_mask(:))/&
        sum(gas_molarconc(:)*vgas_mask(:))

do k=1,num_levels

         Tlay(:,:,num_levels-k+1)  = pt(:,:,k)
         play(:,:,num_levels-k+1)  = pfull(:,:,k) 
!         f_i(:,:,num_levels-k+1,1) = q(:,:,k)
         if(vgas .ne. 0) f_i(:,:,num_levels-k+1,vgas) = q(:,:,k) /&
                (mu_i(vgas)/mu_dry+(1-mu_i(vgas)/mu_dry)*q(:,:,k))
         !from q to molar concentration
enddo
if (vgas==0) then
   mu_avg(:,:,:) = mu_dry
else
   mu_avg(:,:,:)=mu_dry*(1-f_i(:,:,:,vgas)) + mu_i(vgas)*f_i(:,:,:,vgas)
endif
do iGas=1,nGas
   if (iGas .ne. vgas) f_i(:,:,:,iGas) = gas_molarconc(iGas)*(1-f_i(:,:,:,vgas))
enddo
!if (master) write(6,*) mu_dry, mu_avg(1,1,:), &
!        iGas_CO2, f_i(1,1,:,iGas_CO2), &
!        iGas_H2O, f_i(1,1,:,iGas_H2O), &
!        vgas
!calculate f_i(molar concentration) from q and mCO2
!mu_dry = ( 28. + mu_i(2)*mCO2) / (1+mCO2) 
!f_i(:,:,:,1) = f_i(:,:,:,1) / &
!        (mu_i(1)/mu_dry + (1 - mu_i(1)/mu_dry)*f_i(:,:,:,1))
!f_i(:,:,:,2) = mCO2 /(1+mCO2) *(1-f_i(:,:,:,1))
!!for 1D comparison, initialize temperature grids
!-----------------------------------------------------------------------
    !-------- set up collision-induced absorption --------
!    if(iGas_CO2.ne.-1 .and. iGas_H2.ne.-1)then
!       call setup_cia('CO2H2_')
       !call setup_cia('N2_H2_')
       !print*,'Using N2-H2 CIA in place of CO2-H2!'
!    endif
!    if(iGas_CO2.ne.-1 .and. iGas_CH4.ne.-1)then
!       call setup_cia('CO2CH4')
       !call setup_cia('N2_CH4')
       !print*,'Using N2-CH4 CIA in place of CO2-CH4!'
!    endif
    if(iGas_CO2.ne.-1)then
       call setup_cia('CO2CO2')
    endif
    if(iGas_H2O.ne.-1)then
       call setup_cia('H2OH2O')
    endif
    if(iGas_O2.ne.-1)then
        call setup_cia('O2_O2_')
    endif

    !-------- read line data --------
    ! line data is only read if calc_sigma(iGas) = .true.
    if (.not. gray_debug) call read_line_data()

    !-------- set up temperature grid on which cross-sections are defined --------
    call setup_T_k_grid(Tlay(1, 1, :),T_k_grid)

!    if (master) write(6,*) T_k_grid(1,:), T_k_grid(30,:)
!    write(6,*) T_k_grid(num_levels,:) !play(:,:,1), Tlay(:,:,1)

    call define_cosa_quad(cosa,ang_wt)

!define ISR 
!if (tidal_lock) then
!tide-lock, solar_constant
!    cosa0(:,:) =  cos(agrid(:,:,2)) * &
!        cos(agrid(:,:,1) - noon_longitude*deg_to_rad)
!    where (cosa0 .lt. 0.0) cosa0 = 0.0
!    Fstel0(:,:) = solar_constant * cosa0(:,:)
!else
!zonal-symmeric
!    cosa0(:,:) = cos(agrid(:,:,2)) !/ pi
!    Fstel0(:,:) = solar_constant * cosa0(:,:) /pi !zero obliquity
!    Fstel0(:,:) = solar_constant/4.*(1+ 0.25* &
!        (1-3*sin(agrid(:,:,2))**2)) !23.5 obliquity
!endif
!substellar point initialization
kappa_ang = 0.
sublon    = 0.
sublat    = asin(sin(obl*deg_to_rad) *cos(kappa_ang -prec*deg_to_rad))
!zenith angle
cosa0(:,:)= sin(agrid(:,:,2))*sin(sublat) &
            + cos(agrid(:,:,2))*cos(sublat) &
            * cos(agrid(:,:,1) -sublon)
where (cosa0 .lt. 0) cosa0 = 0.0

call orbital_dis(kappa_ang, dist)
Fstel0(:,:) = starlum/(4*pi*dist*dist) *cosa0(:,:)

!gsca = 0.85d0
!call backscatter(gsca, bL_cld)
!do i=1,size(cosa0,1)
!   do j=1,size(cosa0,2)
!      call angular_backscatter(cosa0(i,j), gsca, b0L_cld(i,j))
!   enddo
!enddo
!print *,  b0L_cld(i,1)

    !-------- set up radiative transfer calculations --------
    call setup_shortwave(T_k_grid,play,ps,f_i,mu_i,nu_sw)
!    if (master) write(6,*) play(:,:,1), Tlay(:,:,1)
    call setup_longwave(T_k_grid,play,f_i,mu_i,nu_lw)

    !-------- define lower Tlay for sigma interpolation --------
!    Tcold = Tlay - dTlay


!    if (master) write(6,*) nu0

!define insolation profile------------------------------
!ss  = sin(lat)
!p2 = (1. - 3.*ss*ss)/4.
!solar(:,:) = 0.25*solar_constant*(1.0 + del_sol*p2 + del_sw * ss)
!cos_lat  = cos(lat(1,:))
!   do i = 1, size(t,1)
!      solar(i,:) = 0.25*solar_constant*(1.0 + del_sol*p2 + del_sw * ss)
!      solar(i,:) = solar_constant*cos_lat/pi
!   enddo
!endif

!tide-lock
!      solar(:,:) = solar_constant *cos(agrid(:,:,2)) * &
!      cos(agrid(:,:,1) - noon_longitude*deg_to_rad)
!      where (solar .lt. 0.0) solar = 0.0

!------------ initialize diagnostic fields ---------------


id_tau_lw_inf  = register_diag_field(mod_name, 'tau_lw_inf',        &
     axes((/1,2,5/)), Time, 'longwave optical thickness','none') 

id_flux_lw   = register_diag_field(mod_name, 'flux_lw',        &
     axes((/1,2,4/)), Time, 'longwave flux(positive upward)','W/m**2') 

id_flux_sw   = register_diag_field(mod_name, 'flux_sw',        &
     axes((/1,2,4/)), Time, 'shortwave flux(positive downward)','W/m**2')

id_sigma_lw   = register_diag_field(mod_name, 'sigma_lw',        &
     axes((/3,5/)), Time, 'longwave abs','') 

id_sigma_di   = register_diag_field(mod_name, 'sigma_di',        &
     axes((/1,2,5/)), Time, 'longwave abs','') 

id_sigma_vi   = register_diag_field(mod_name, 'sigma_vi',        &
     axes((/1,2,5/)), Time, 'longwave abs','') 

id_OLRnu = register_diag_field(mod_name, 'OLRnu',        &
     axes((/1,2,5/)), Time, 'OLR spectrum','W/m**2/cm-1') 

id_OLR   = register_diag_field(mod_name, 'OLR',        &
     axes((/1,2/)), Time, 'OLR','W/m**2') 

id_cosa  = register_diag_field(mod_name, 'cosa', &
     axes((/1,2/)), Time, 'cos of incident angle','none')

id_solar = register_diag_field(mod_name, 'solar', &
     axes((/1,2/)), Time, 'insolation','W/m**2')

id_solarnu = register_diag_field(mod_name, 'solarnu', &
     axes((/1,2,5/)), Time, 'insolation spectrum','W/m**2/cm-1')

id_OSRnu = register_diag_field(mod_name, 'OSRnu',        &
     axes((/1,2,5/)), Time, 'OSR spectrum','W/m**2/cm-1') 

id_OSR   = register_diag_field(mod_name, 'OSR',        &
     axes((/1,2/)), Time, 'OSR','W/m**2') 

id_hrlw  = register_diag_field(mod_name, 'hrlw',        &
     axes((/1,2,3/)), Time, 'longwave heating rate','K/sec') 

id_hrsw  = register_diag_field(mod_name, 'hrsw',        &
     axes((/1,2,3/)), Time, 'shortwave heating rate','K/sec') 

id_Fdn_sw   = register_diag_field(mod_name, 'Fdn_sw',        &
     axes((/1,2,4/)), Time, 'shortwave flux(downward diffusive)','W/m**2')

id_Fup_sw   = register_diag_field(mod_name, 'Fup_sw',        &
     axes((/1,2,4/)), Time, 'shortwave flux(upward diffusive)','W/m**2')

id_Fdir_sw   = register_diag_field(mod_name, 'Fdir_sw',        &
     axes((/1,2,4/)), Time, 'shortwave flux(direct beam)','W/m**2')

!if (id_cosa>0)     used = send_data(id_cosa , cosa0(:,:), Time)
!if (id_solar>0)    used = send_data(id_solar, Fstel0(:,:), Time)
!if (id_solarnu>0)  used = send_data(id_solarnu, F_stel_dn(:,:,:), Time)
if (id_sigma_lw>0) used = send_data(id_sigma_lw, sigma_lw_ar(:,:,2,1) , Time) !npz, nS, nGas, nTem

return
end subroutine radiance_init


! ==================================================================================

subroutine radiance_update (is, ie, js, je, num_levels, &
        Time_diag, agrid, q, q_cld, pt, ts, bucket, phalf, pfull, dTrad, dTs_rad )

integer, intent(in)                 :: is, ie, js, je, num_levels
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:,:)  :: agrid
real, intent(in), dimension(:,:,:)  :: q, pt, pfull, q_cld
real, intent(in), dimension(:,:,:)  :: phalf
real, intent(in), dimension(:,:)    :: ts, bucket
real, intent(out), dimension(:,:,:) :: dTrad
real, intent(out), dimension(:,:)   :: dTs_rad

integer :: iGas, iAng, iLev, iLay, nLay, i, j

!master = (gid == masterproc)

nLay = num_levels
! reverse the order, 1st is surface now

do iLay=1,nLay
         Tlay(:,:,nLay-iLay+1)  = pt(:,:,iLay)
         play(:,:,nLay-iLay+1)  = pfull(:,:,iLay) 
!         f_i(:,:,nLay-iLay+1,1) = q(:,:,iLay)
         plev(:,:,nLay-iLay+2)  = phalf(:,:,iLay) 

         if(vgas .ne. 0) f_i(:,:,nLay-iLay+1,vgas) = q(:,:,iLay) /&
                (mu_i(vgas)/mu_dry+(1-mu_i(vgas)/mu_dry)*q(:,:,iLay))

         f_cloud(:,:,nLay-iLay+1)= q_cld(:,:,iLay)
         !lnplev(i-is+1,j-js+1,nLay-iLay+2)= peln(i,iLay,j) 
enddo
         plev(:,:,1) = phalf(:,:,nLay+1) !peln(i,nLay+1,j)
         t_surf(:,:)  = ts(:,:)

!compute Tlev
dp(:,:,:) = plev(:,:,1:nLay) - plev(:,:,2:nLay+1)  
lnplev = log(plev)
lnplay = log(play)
do i=1,size(Tlay,1)
   do j=1,size(Tlay,2)
      call interp(nLay, lnplay(i,j,:), Tlay(i,j,:), &
                nLay+1, lnplev(i,j,:), Tlev(i,j,:))
      Tlev(i,j,1) = 2.*Tlay(i,j,1) - Tlev(i,j,2)
      Tlev(i,j,nLay+1) = 2.*Tlay(i,j,nLay) - Tlev(i,j,nLay)
   enddo
enddo

if (vgas==0) then
   mu_avg(:,:,:) = mu_dry
else
   mu_avg(:,:,:)=mu_dry*(1-f_i(:,:,:,vgas)) + mu_i(vgas)*f_i(:,:,:,vgas)
endif
do iGas=1,nGas
   if (iGas .ne. vgas) f_i(:,:,:,iGas) = gas_molarconc(iGas)*(1-f_i(:,:,:,vgas))
enddo
! if (master) write(6,*) mu_dry, mu_avg(1,1,:), &
!         iGas_CO2, mu_i(iGas_CO2),  &
!         iGas_H2O, mu_i(iGas_H2O),  &
!         iGas_N2 , mu_i(iGas_N2 ), f_i(:,:,1,vgas), &
!         vgas
!!for 1D comparison
!t_surf(:,:) = 281.
!do iLay=1,nLay
!   Tlay(:,:,iLay) = t_surf(:,:)* (play(:,:,iLay)/plev(:,:,1))**kappa
!   Tlev(:,:,iLay+1) = t_surf(:,:)* &
!        (plev(:,:,iLay+1)/plev(:,:,1))**kappa
!enddo
!Tlev(:,:,1)  = t_surf(:,:)
!do i=1,size(Tlay,1)
!   do j=1,size(Tlay,2)
!      f_i(i,j,:,1) = (/1.63520966e-02,   1.57950129e-02,   1.45766120e-02, &
!         1.17701637e-02,   9.78253968e-03,   7.57087674e-03, &
!         5.25058620e-03,   3.12429061e-03,   1.64496037e-03, &
!         7.98964349e-04,   3.72633222e-04,   1.57618095e-04, &
!         5.96557402e-05,   2.05953020e-05,   6.54622499e-06, &
!         2.10978124e-06,   6.60884609e-07,   2.11125069e-07, &
!         6.09305033e-08,   1.42561145e-08,   1.28652600e-09, &
!         2.74697626e-10,   1.89126770e-10,   1.23416791e-10, &
!         6.35112310e-11,   1.10281111e-10/) 
!   enddo
!enddo
!f_i(:,:,:,2) = mCO2
!mu_avg(:,:,:)= 28.
!substellar point update
call find_substellar(kappa_ang, sublon, kappa_new, sublon_new, sublat)
kappa_ang = kappa_new *1.
sublon    = sublon_new *1. 
!zenith angle
cosa0(:,:)= sin(agrid(:,:,2))*sin(sublat) &
            + cos(agrid(:,:,2))*cos(sublat) &
            * cos(agrid(:,:,1) -sublon)
where (cosa0 .lt. 0) cosa0 = 0.0

call orbital_dis(kappa_ang, dist)
Fstel0(:,:) = starlum/(4*pi*dist*dist) *cosa0(:,:)

    call init_stellar_spectrum(sigma_sw)
    do i=1,size(play,1)
       do j=1,size(play,2)
          if (cosa0(i,j) .gt. 0) then
             F_stel_dn(i,j,:) = sigma_sw(:)*Fstel0(i,j)/sum(sigma_sw(:)*dnu_sw) ! stellar spectral flux [W/m2/cm^-1]
             I_stel_dn(i,j,:) = F_stel_dn(i,j,:)/(Omega_stel*cosa0(i,j))     ! stellar spectral irradiance [W/m2/cm^-1/sr]
    ! from classic flux defn. (e.g. Goody & Yung)
          else
             F_stel_dn(i,j,:) = 0.
             I_stel_dn(i,j,:) = 0.
          endif
       enddo
    enddo

    !-------- temperature grid interpolation calculation --------
    T_lin_weight = -1.0e8
    if(nTem>1) call calculate_T_k_grid(Tlay(:,:,:),iT1,T_lin_weight)

    !-------- visible radiative transfer calculation --------
    
    OSRnu(:,:,:)        = 0.0d0
    OSR(:,:)            = 0.0d0
    GSR(:,:)            = 0.0d0
    Flev_up_sw(:,:,:)   = 0.0d0
    Flev_dn_sw(:,:,:)   = 0.0d0
    Flev_dir_sw(:,:,:)  = 0.0d0

    !if(verbose)then
    !   write(*,*) ' Calculating shortwave direct beam attenuation.'
    !endif
    sw_angle_loop: do iAng = 1, nAng ! integration over propagation angle

!       if(verbose)then
!          write(*,'(a44,f8.2,a9)') ' Calculating shortwave with emission angle = ', &
!               acos(cosa(iAng))*180.0d0/pi, ' degrees.'
!       endif
       call calculate_shortwave(iT1,cosa(iAng),T_lin_weight,t_surf,dp,f_i,&
               f_cloud,mu_avg,kappa_gray, bucket, &
               OSRnu_temp,OSR_temp,GSR_temp,Ilev_up_sw,Ilev_dn_sw,Flev_dir_sw)

       ! c.f. e.g. Stamnes et al. (2000), DISORT tech. report, eqn. (9a)
       OSR        = OSR        + 2*pi*OSR_temp  *cosa(iAng)*ang_wt(iAng)
       OSRnu      = OSRnu      + 2*pi*OSRnu_temp*cosa(iAng)*ang_wt(iAng)
       GSR        = GSR        + 2*pi*GSR_temp  *cosa(iAng)*ang_wt(iAng)
       Flev_up_sw = Flev_up_sw + 2*pi*Ilev_up_sw*cosa(iAng)*ang_wt(iAng)
       Flev_dn_sw = Flev_dn_sw + 2*pi*Ilev_dn_sw*cosa(iAng)*ang_wt(iAng)

    end do sw_angle_loop

    GSR = GSR + Flev_dir_sw(:,:,1)
    !GSRnu      = GSRnu_temp
    ! absorbed stellar radiation = ISR - OSR [W/m2]
    !ASR = ISR - OSR

    !-------- IR radiative transfer calculation --------
    
    OLRnu(:,:,:)        = 0.0d0
    OLR(:,:)            = 0.0d0
    Flev_up_lw(:,:,:)   = 0.0d0
    Flev_dn_lw(:,:,:)   = 0.0d0

    lw_angle_loop: do iAng = 1, nAng ! integration over propagation angle

       !if(verbose)then
       !   write(*,'(a44,f8.2,a9)') ' Calculating longwave with emission angle = ', &
       !        acos(cosa(iAng))*180.0d0/pi, ' degrees.'
       !endif

       call calculate_longwave(iT1,cosa(iAng),Tlay,Tlev,t_surf,T_lin_weight,play,dp,f_i,&
               f_cloud,mu_avg, &
               kappa_gray,OLRnu_temp,OLR_temp,Ilev_up_lw,Ilev_dn_lw)

       OLR        = OLR        + 2*pi*OLR_temp  *cosa(iAng)*ang_wt(iAng)
       OLRnu      = OLRnu      + 2*pi*OLRnu_temp*cosa(iAng)*ang_wt(iAng)
       Flev_up_lw = Flev_up_lw + 2*pi*Ilev_up_lw*cosa(iAng)*ang_wt(iAng)
       Flev_dn_lw = Flev_dn_lw + 2*pi*Ilev_dn_lw*cosa(iAng)*ang_wt(iAng)

    end do lw_angle_loop

    !-------- heating rate calculation --------
    do ilay=1,nLay
       iLev = iLay + 1
       dF_lw(:,:,iLay) = (Flev_up_lw(:,:,iLev) - Flev_up_lw(:,:,iLev-1)) - &
       (Flev_dn_lw(:,:,iLev) - Flev_dn_lw(:,:,iLev-1))
       dF_sw(:,:,iLay) = (Flev_up_sw(:,:,iLev) - Flev_up_sw(:,:,iLev-1)) - &
       (Flev_dn_sw(:,:,iLev) - Flev_dn_sw(:,:,iLev-1)) - &
       (Flev_dir_sw(:,:,iLev) - Flev_dir_sw(:,:,iLev-1))
    end do
    dTdt_lw(:,:,:) = -(grav/cp_air )*dF_lw(:,:,:)/dp(:,:,:)
    dTdt_sw(:,:,:) = -(grav/cp_air )*dF_sw(:,:,:)/dp(:,:,:)
    dTdt(:,:,:)    = dTdt_lw(:,:,:) + dTdt_sw(:,:,:)
!    dTdt(:,:,:)    = -dF_lw(:,:,:)
    !net_surf_sw_down(:,:)       = Flev_dn_sw(:,:,1) !solar
    !surf_lw_down(:,:)           = Flev_dn_lw(:,:,1)

    dt_surf_rad =  Flev_dn_lw(:,:,1) - Flev_up_lw(:,:,1) + &
        Flev_dn_sw(:,:,1) - Flev_up_sw(:,:,1) + &
        Flev_dir_sw(:,:,1)
               !    dt_surf_rad = (net_surf_sw_down + Flev_dn_lw(:,:,1) - &
!                Flev_up_lw(:,:,1)) / rho_cp
!    dt_surf_rad = (surf_lw_down + net_surf_sw_down  &
!                     - 5.67e-8*t_surf **4) &
!                     / rho_cp !/ mld !eff_heat_capacity
!update temperature
!Tlay   = Tlay + dTdt *dt_atmos
!t_surf = t_surf + dt_surf_rad *dt_atmos 

!reverse the order again

do iLay=1,nLay
   dTrad(:,:,iLay) = dTdt(:,:,nLay-iLay+1)  
enddo
   dTs_rad(:,:) = dt_surf_rad(:,:) 

    !-------- report results --------
    !if(report_results)then
    !   write(*,*) '---------------------------------'
    !   write(*,*) 'Results:'
    !   write(*,'(a9,f14.7,a5)') '  ISR  = ',ISR,' W/m2'
    !   write(*,'(a9,f14.7,a5)') '  OLR  = ',OLR,' W/m2'
    !   write(*,'(a9,f14.7,a5)') '  ASR  = ',ASR,' W/m2'
    !   write(*,'(a9,f14.7)')    '  OLRe = ',OLR/(stefan*Ts**4)
    !   write(*,'(a9,f14.7)')    '  Apla = ',OSR/ISR
    !   write(*,'(a9,f14.7)')    '  Areq = ',1.0d0 - OLR/ISR
       ! Areq is the albedo required for thermal equilibrium
    !endif

!------- downward lw flux surface -------
!      if ( id_lwdn_sfc > 0 ) then
!          used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
!      endif

!OLR = Tlev(:,:,1)
if (id_hrsw>0) used=send_data(id_hrsw, ( dtau_sw_ray(:,:,:,188) +dtau_sw_cld(:,:,:,188)) /&
            dtau_sw_i(:,:,:,188) , Time_diag)
!if (id_hrsw>0) used=send_data(id_hrsw, dTdt_sw, Time_diag)
!if (id_hrsw>0) used=send_data(id_hrsw, dTdt_sw, Time_diag)
if (id_hrlw>0) used=send_data(id_hrlw, dTdt_lw, Time_diag)

if (id_tau_lw_inf>0)    used = send_data(id_tau_lw_inf, &
        tau_lw_inf, Time_diag)
if (id_flux_lw>0)       used = send_data(id_flux_lw, &
        Flev_up_lw-Flev_dn_lw, Time_diag)
if (id_flux_sw>0)       used = send_data(id_flux_sw, &
        Flev_dir_sw+Flev_dn_sw-Flev_up_sw, Time_diag)
!if (id_tau_lw_inf>0)    used = send_data(id_tau_lw_inf, &
!        iT1, Time_diag)

!if (master) write(6,*) Tlay(1,1,:), Tlev(1,1,:)
!if (master) write(6,*) t_surf(:,:)
if (id_sigma_di>0) used = send_data(id_sigma_di, &
!        sigma_d_i(:,:,nLay,:) , Time_diag)
        sigma_d_i(:,:,1,:) , Time_diag)
if (id_sigma_vi>0) used = send_data(id_sigma_vi, &
        sigma_v_i(:,:,1,:) , Time_diag)

if (id_OLRnu>0) used = send_data(id_OLRnu, OLRnu, Time_diag)
if (id_OLR>0)   used = send_data(id_OLR, OLR, Time_diag)
!if (id_OLR>0)   used = send_data(id_OLR, t_surf(:,:), Time_diag)
if (id_OSRnu>0) used = send_data(id_OSRnu, OSRnu, Time_diag)
if (id_OSR>0)   used = send_data(id_OSR, OSR, Time_diag)

if (id_Fdn_sw>0) used=send_data(id_Fdn_sw, Flev_dn_sw, Time_diag)
if (id_Fup_sw>0) used=send_data(id_Fup_sw, Flev_up_sw, Time_diag)

if (id_Fdir_sw>0) used=send_data(id_Fdir_sw, Flev_dir_sw, Time_diag)

if (id_cosa>0)     used = send_data(id_cosa , cosa0(:,:), Time_diag)
if (id_solar>0)    used = send_data(id_solar, Fstel0(:,:), Time_diag)
if (id_solarnu>0)  used = send_data(id_solarnu, F_stel_dn(:,:,:), Time_diag)

!call get_CIA_cross_section(iGas_CO2,nu_lw,sigma_CIA,Tlay(1,1,1),play(1,1,1),f_i(1,1,1,:))
!if (master) write(6,*) sigma_CIA(:)

return
end subroutine radiance_update

! ==================================================================================

                                                                      
subroutine radiance_end
                                                                                                      
!deallocate (b, tdt_rad, tdt_sw, entrop_rad) 


end subroutine radiance_end

! ==================================================================================


 subroutine read_line_data()

    ! read all line data
    integer :: kk, il, ierr, iGas

    integer :: mol_temp, iso_temp      ! molecule number
    real ::  nu0_temp, Sref_temp, einstein_temp,gamair_temp, gamself_temp,Egnd_temp,ncoeff_temp,  delta_temp

    ! read lines
    if (master) write(6,*) "creating line data..."

    molec_loop : do iGas = 1, nGas
!    if(calc_sigma(iGas))then
       if(iGas==iGas_H2O)then
          if(HITEMP)then
             open(unit=111,file=trim(datadir)//'01_HITEMP2010_red2.par')
          else
             open(unit=111,file=trim(datadir)//'01_hit12.par')
          endif
       elseif(iGas==iGas_CO2)then
          if(HITEMP)then
             open(unit=111,file=trim(datadir)//'02_HITEMP2010_red2.par')
          else
             open(unit=111,file=trim(datadir)//'02_hit12.par')
          endif
!       elseif(iGas==iGas_O3)then
!          open(unit=111,file=trim(datadir)//'03_hit12.par')
!       elseif(iGas==iGas_CH4)then
!          open(unit=111,file=trim(datadir)//'06_hit12.par')
!       elseif(iGas==iGas_SO2)then
!          open(unit=111,file=trim(datadir)//'09_hit12.par')
!       elseif(iGas==iGas_H2S)then
!          open(unit=111,file=trim(datadir)//'31_hit12.par')
       else
          !write(6,*) 'Note: spectral dataset not found for iGas = ',iGas,' .'
          cycle molec_loop
       end if

       il = 1
       read_loop : do 
          if(HITEMP)then
             read(111,'(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2)',iostat=kk) mol_temp, iso_temp, &
                  nu0_temp, Sref_temp, einstein_temp, gamair_temp, gamself_temp, Egnd_temp, ncoeff_temp!, delta_temp
          else
             read(111,'(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2, f8.2)',iostat=kk) mol_temp, iso_temp, &
                  nu0_temp, Sref_temp, einstein_temp, gamair_temp, gamself_temp, Egnd_temp, ncoeff_temp, delta_temp
          endif

          !if(Sref_temp>1.0d-25)then ! don't allow weak lines
          if(Sref_temp>1.0d-25)then ! don't allow weak lines
          !if(Sref_temp>1.0d-30)then ! don't allow weak lines
          !if(Sref_temp>1.0d-35)then ! don't allow weak lines
             mol(iGas,il)      = mol_temp
             iso(iGas,il)      = iso_temp
             nu0(iGas,il)      = nu0_temp
             Sref(iGas,il)     = Sref_temp
             einstein(iGas,il) = einstein_temp
             gam_air(iGas,il)  = gamair_temp
             gam_self(iGas,il) = gamself_temp
             Egnd(iGas,il)     = Egnd_temp
             ncoeff(iGas,il)   = ncoeff_temp
             delta(iGas,il)    = delta_temp
             il = il + 1
          endif

          
          if(kk == -1)then
             exit read_loop
          endif
          if(il>nlinesMAX)then
             exit read_loop
          endif

       end do read_loop
       close(111)
       
       nlines(iGas) = il - 1

       ! special extra to add unassigned near-IR methane lines from
       ! Beguier+, JQRST (2015)
!       if(use_beguier_2015 .and. iGas==iGas_CH4)then
!          write(*,*) 'Adding Beguier+ (2015) unassigned CH4 lines in near-IR...'
!          open(unit=112,file=trim(datadir)//'beguier_2015.par')
!          do il=nlines(iGas)+1,nlines(iGas)+nlines_beguier
!             read(112,*) nu0_temp, Sref_temp
!             mol(iGas,il)      = 6
!             iso(iGas,il)      = 1
!             nu0(iGas,il)      = nu0_temp
!             Sref(iGas,il)     = Sref_temp
!             einstein(iGas,il) = -1.0e10
!             gam_air(iGas,il)  = 1.8d-2 ! from Beguier+ (2015) Section 3.1.
!             gam_self(iGas,il) = 0.0d0
!             Egnd(iGas,il)     = 1.0d0
!             ncoeff(iGas,il)   = 0.0d0 ! TO BE DETERMINED
!             delta(iGas,il)    = 0.0d0
!          end do
!          close(112)

!          nlines(iGas) = il - 1
!       end if


!    end if
    end do molec_loop

  end subroutine read_line_data

  ! ==================================================================================

  subroutine line_strength_width(mu_i,p,p_i,T,iGas)

    ! calculate line strengths, widths and frequency pressure shifts
    ! for a given species at a given temperature

    real, intent(in) :: T    ! temperature [K]
    real, intent(in) :: p    ! total pressure [Pa]
    real, intent(in) :: p_i  ! partial pressure [Pa]
    real, intent(in) :: mu_i ! molar mass [g/mol]
    integer, intent(in) :: iGas

    real :: Tfact, Tfact1       ! temperature-dependent factors []
    real :: QrefQ               ! total internal partition sum
    !real(8) gi                  ! state independent degeneracy factor
    
    integer :: il

    logical, parameter :: do_air_shift = .false. ! calculate pressure-shifted line center
    logical, parameter :: verbose      = .false. ! print more info for diagnostic

    if(T<20.0d0)then
       !print *, T, p, p_i
       call error_mesg('radiance:','In line_strength_width T = exiting',FATAL)
    end if

    if(iGas==iGas_O2 .or. iGas==iGas_N2)then ! todo: generalize this
!       if(verbose) write(*,*) 'Treating ', gas_name(iGas), ' as vib-rot inactive.'
    else

!       if(verbose) write(*,*) "calculating line strengths for ",gas_name(iGas),"..."

       
       read_loop : do il = 1, nlines(iGas)
          
          ! calculate (air) pressure-shifted nu
!          if(do_air_shift)then
             !nu_shift(iGas,il) = nu0(iGas,il) + delta(il)*p_i(iGas)/atm
!          end if
          
          ! calculate Lorentz and Doppler linewidths [cm^-1]
          ! Lorentz linewidth is the Lorentz HWHM
          ! Doppler linewidth is the Doppler HWHM / sqrt(ln(2))
          gam_lore(iGas,il) = (Tref/T)**ncoeff(iGas,il) &
                               *(gam_air(iGas,il)*(p - p_i)/atm + gam_self(iGas,il)*p_i/atm)
          gam_dopp(iGas,il) = nu0(iGas,il)*sqrt(2*Rstar*T/(mu_i/1.0d3))/c
          ! note mu_i not mu_avg here because it is the mean speed of the molecule doing the
          ! absorbing that counts.
          
          ! get Total Internal Partition Sum (TIPS) ratio Qref / Q
          ! simple approach based on BD_TIPS data
          ! robust up to about 1000 K for H2O and CO2
          ! up to 500-600 K for CH4
          ! up to about 350 K for O3
          if(iGas==iGas_H2O .or. iGas==iGas_CO2 .or. iGas==iGas_O3 .or. iGas==iGas_CH4)then
             QrefQ = (Tref/T)**1.5d0
          else
             call error_mesg('radiance: ','Molecule Q(T) needs to be assessed still!',FATAL)
          end if

          
          ! get actual line intensity
          Tfact1          = (1.0d0 - exp(-c2*nu0(iGas,il)/T)) &
                              / (1.0d0 - exp(-c2*nu0(iGas,il)/Tref))
          Tfact           = exp( -c2*Egnd(iGas,il)*(1.0d0/T - 1.0d0/Tref) )
          Strue(iGas,il)  = Sref(iGas,il)*Tfact*Tfact1*QrefQ

          ! no temperature scaling for Beguier+ (2015) lines
          if(use_beguier_2015 .and. iGas.eq.iGas_CH4 .and. il>(nlines(iGas)-nlines_beguier)) Strue(iGas,il) = Sref(iGas,il)
          
       end do read_loop
       
       nlines(iGas) = il - 1

    end if
    
  end subroutine line_strength_width

 
  ! ==================================================================================


  subroutine get_line_abs_cross_section(iGas,nu,sigma,T)
    
    ! calculate line absorption cross-section vs. nu for a given species
    ! in a given wavenumber range

    integer, intent(in)  :: iGas
    real, intent(in)  :: nu(nS)     ! wavenumber [cm^-1]
    real, intent(in)  :: T          ! temperature [K]: for CO2 chi-factor calculation
    real, intent(out) :: sigma(nS)  ! total absorption cross-section [cm2/molec]

    logical :: mask(nS)
    integer :: iS, il
    real :: f_temp, gam_temp
    
    real :: nu1  ! start wavenumber [cm^-1]
    real :: nu2  ! finish wavenumber [cm^-1]
    
    logical, parameter :: use_voigt   = .true.  ! use voigt function for line profiles
    logical, parameter :: debug       = .false. ! print additional text for diagnostic

!    if(debug) write(*,*) "calculating absorption cross sections for ",gas_name(iGas),"..."

    sigma(:) = 0.0d0

    nu1 = nu(1)
    nu2 = nu(nS)

    if(iGas==iGas_O2 .or. iGas==iGas_N2)then ! todo: transfer update from generate_kmatrix

    else
       ! add the lines, one by one
       ! for now we use all isotopes (with Earth-like abundances!)
       line_loop : do il = 1, nlines(iGas)

          ! this is the very slow part
          trunc : if(nu0(iGas,il)>nu1-deltanu_trunc(iGas) .and. nu0(iGas,il)<nu2+deltanu_trunc(iGas))then

             mask = abs(nu - nu0(iGas,il)) < deltanu_trunc(iGas)

             write_sig_loop : do iS = 1, nS

                ! update sigma only when mask = T
                if(mask(iS))then

                   if(iGas==iGas_CO2)then
                      gam_temp = gam_lore(iGas,il) *chi_factor(nu(iS),nu0(iGas,il),T)
                      ! check but I believe this is correct
                   else
                      gam_temp = gam_lore(iGas,il) ! no sublorentzian lineshift
                   end if

                   if(use_voigt)then
                      call voigt(gam_temp, gam_dopp(iGas,il), nu(iS)-nu0(iGas,il), f_temp)
                   else
                      call lorentz(nu(iS)-nu0(iGas,il),gam_temp,f_temp)
                      !call doppler(nu(iS)-nu0(iGas,il),gam_dopp(iGas,il),f_temp)
                   endif

                   if(debug)then
                      call voigt(0.0d0,1.0d0,0.0d0,f_temp)
                      ! f_V(0) when Lorentz HWHM = 0: should yield f_D(0) = 1/sqrt(pi)
                      print*,f_temp
                      call voigt(1.0d0,1.0d-8,0.0d0,f_temp)
                      ! f_V(0) when Doppler HWHM -> 0: should yield f_L(0) = 1/pi
                      print*,f_temp
                      stop
                   endif
                   
                   sigma(iS) = sigma(iS) + Strue(iGas,il)*f_temp
                endif
             end do write_sig_loop
          end if trunc
       end do line_loop

    end if
    
  end subroutine get_line_abs_cross_section

  ! ==================================================================================
  ! ==================================================================================

  subroutine get_CIA_cross_section(iGas,nu,sigma_CIA,T,p,f_i)

    ! calculate cross-section vs. nu for a given species
    ! in a given wavenumber range

    implicit none

    integer, intent(in)  :: iGas
    real(8), intent(in)  :: nu(nS)        ! wavenumber [cm^-1]
    real(8), intent(in)  :: T             ! temperature [K]
    real(8), intent(in)  :: p             ! pressure [Pa]
    real(8), intent(in)  :: f_i(nGas)     ! molar concentration [mol/mol]
    real(8), intent(out) :: sigma_CIA(nS) ! CIA absorption cross-section [cm2/molecule of air]

    integer iS
    real(8) sigma_temp                  ! absorption cross-section [cm2/molecule of air]

    sigma_CIA(:) = 0.0d0

    !--------------------------------------------------------------
    ! add any relevant collision-induced absorption (longwave only)

    if(iGas==iGas_H2O)then
       do iS = 1, nS
          sigma_temp = 0.0d0
          call interpolateH2Ocont_CKD(nu(iS),T,p*f_i(iGas_H2O),p*(1.0d0-f_i(iGas_H2O)),sigma_temp,.false.)
          sigma_CIA(iS) = sigma_temp
       end do
    end if

    if(iGas==iGas_CO2)then
       do iS = 1, nS
          sigma_temp = 0.0d0
          if(nu(iS)>0.0d0 .and. nu(iS)<500.0d0)then
             call calculate_CO2CO2_cia(1,T,nu(iS),p*f_i(iGas_CO2)/atm,sigma_temp) ! induced dipole             
          elseif(nu(iS)>700.0d0 .and. nu(iS)<2000.0d0) then
             call calculate_CO2CO2_cia(2,T,nu(iS),p*f_i(iGas_CO2)/atm,sigma_temp) ! dimer
          end if
          sigma_CIA(iS) = sigma_temp * f_i(iGas_CO2) !different from O2 due to nomorlization
       end do
    end if

     if(iGas==iGas_O2)then
        do iS = 1, nS
           sigma_temp = 0.0d0
           if(nu(iS)>1150.0d0 .and. nu(iS)<1950.0d0)then
              call calculate_O2O2_cia(T, nu(iS),p*f_i(iGas_O2)/atm,f_i(iGas_O2),sigma_temp) 
           elseif(nu(iS)>7450.0d0 .and. nu(iS)<8487.0d0) then
              call calculate_O2O2_cia(T, nu(iS),p*f_i(iGas_O2)/atm,f_i(iGas_O2),sigma_temp) 
           elseif(nu(iS)>9001.0d0 .and. nu(iS)<9997.0d0) then
              call calculate_O2O2_cia(T, nu(iS),p*f_i(iGas_O2)/atm,f_i(iGas_O2),sigma_temp) 
           elseif(nu(iS)>12600.0d0 .and. nu(iS)<13839.0d0) then
              call calculate_O2O2_cia(T, nu(iS),p*f_i(iGas_O2)/atm,f_i(iGas_O2),sigma_temp) 
           elseif(nu(iS)>14996.0d0 .and. nu(iS)<29790.0d0) then
              call calculate_O2O2_cia(T, nu(iS),p*f_i(iGas_O2)/atm,f_i(iGas_O2),sigma_temp) 
           end if
           sigma_CIA(iS) = sigma_temp
        end do
     end if

!    if(iGas==iGas_H2 .and. iGas_CO2.ne.-1)then ! H2-CO2 CIA gets called once, for H2
!       do iS = 1, nS
!          sigma_temp = 0.0d0
!          if(nu(iS)>10.0d0 .and. nu(iS)<1880.0d0)then
!             !call calculate_cia('N2_H2_',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_H2),sigma_temp,.false.)
!             !print*,'using N2-H2 for CO2-H2!'
!             call calculate_cia('CO2H2_',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_H2),sigma_temp,.false.)
!          end if
!          sigma_CIA(iS) = sigma_temp
!       end do
!    end if

    if(iGas==iGas_CH4 .and. iGas_CO2.ne.-1)then ! CO2-CH4 gets called once, for CH4

       do iS = 1, nS
          sigma_temp = 0.0d0
          if(nu(iS)>10.0d0 .and. nu(iS)<1370.0d0)then
             !print*,'ignoring CO2-CH4 CIA!'
             !call calculate_cia('N2_CH4',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_CH4),sigma_temp,.false.)
             call calculate_cia('CO2CH4',nu(iS),T,p,p*f_i(iGas_CO2),p*f_i(iGas_CH4),sigma_temp,.false.)
          end if
          sigma_CIA(iS) = sigma_temp
       end do
    end if

  end subroutine get_CIA_cross_section

  ! ==================================================================================
  real(8) function chi_factor(nu_temp,nu_L,T)

    ! calculate sublorentzian line profile chi factor using Perrin & Hartman (1989) assumptions

    ! to be double-checked

    implicit none

    real(8), intent(in) :: T       ! temperature [K]
    real(8), intent(in) :: nu_L    ! line center wavenumber [cm^-1]
    real(8), intent(in) :: nu_temp ! wavenumber [cm^-1]

    ! parameters for the empirical scheme
    real(8), parameter, dimension(3) :: alpha = [0.0888d0,  0.0d0,     0.0232d0]
    real(8), parameter, dimension(3) :: beta  = [-0.160d0,  0.0526d0,  0.0d0]
    real(8), parameter, dimension(3) :: epsi  = [0.00410d0, 0.00152d0, 0.0d0]
    real(8), parameter, dimension(3) :: sigPH = [3.0,       30.0,      120.0] ! cm^-1

    real(8), dimension(3) :: B
    real(8) deltanu ! wavenumber separation [cm^-1]
    
    B          = alpha + beta*exp(-epsi*T)
    chi_factor = 1.0d0;

    deltanu   = abs(nu_temp - nu_L)
    if(deltanu<sigPH(1))then
       chi_factor = 1.0d0
    elseif(deltanu>sigPH(1) .and. deltanu<sigPH(2))then
       chi_factor = exp(-B(1)*(deltanu - sigPH(1)))
    elseif(deltanu>sigPH(2) .and. deltanu<sigPH(3))then
       chi_factor = exp(-B(1)*(sigPH(2) - sigPH(1))-B(2)*(deltanu - sigPH(2)))
    else
       chi_factor = exp(-B(1)*(sigPH(2) - sigPH(1)) - B(2)*(sigPH(3) - sigPH(2)) &
            - B(3)*(deltanu - sigPH(3)))
    end if

    return
  end function chi_factor

  ! ==================================================================================
  

  subroutine setup_longwave(T_k_grid,play,f_i,mu_i,nu_lw_out)

    real, intent(in)  :: T_k_grid(:,:)        ! temperature in layer  [K]
    real, intent(in)  :: play(:,:,:)        ! pressure layers [Pa]
    real, intent(in)  :: f_i(:,:,:,:)    ! species molar concentration [mol/mol]
    real, intent(in)  :: mu_i(nGas)        ! molar mass of species [g/mol]
    real, intent(out) :: nu_lw_out(nS)     ! longwave wavenumber [cm^-1]

    real :: nu_lw_check                       ! longwave wavenumber [cm^-1]
    integer ::  iGas, iLay, iS, iTem
    !------- subroutine options --------------  
    logical, parameter   :: verbose = .true.

 !   ! read namelist 
 !   open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
 !   if (ierr/=0) then
 !      print*, 'Cannot find required input.nml file, aborting.'
 !      call abort
 !   else
 !      read(10,longwave_nml)
 !      close(10)
 !   endif
    
    ! create fine spectral grid and initialize abs array
    nu_lw(1)  = nu_lw1
    dnu_lw = (nu_lw2-nu_lw1)/dble(nS-1)  ! Wavenumber interval [cm^-1]
    do iS = 2, nS
       nu_lw(iS) = nu_lw(iS-1) + dnu_lw
    end do
    nu_lw_out = nu_lw

    ! calculate absorption cross section at each level

    sigma_lw(:)                = 0.0d0
    tau_lw_inf(:,:,:)          = 0.0d0
    sigma_lw_ar(:,:,:,:)       = 0.0d0
    sigma_total_dry(:,:,:)     = 0.0d0
    sigma_t_i(:,:,:,:)         = 0.0d0
    dtau_lw(:,:,:,:)           = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low

    gray_debug_if : if(gray_debug)then

       ! don't calculate anything in this case

    else

       big_species_loop : do iGas = 1, nGas
          
          ! read gas absorption cross section data from file if requested
!          if(.not.calc_sigma(iGas))then

             ! no need to check the number of temperatures for each layer match
             ! this was done in setup_shortwave
             
!             ! check the spectral arrays match
!             open(unit=2,file='saved_sigma_data/nu_lw.dat')
!             do iS = 1, nS
!                read(2,*) nu_lw_check
!                if(abs(nu_lw_check-nu_lw(iS))>1.0d-8)then
!                   write(*,*) 'Error: nu_lw data in saved_sigma_data does not match!'
!                   print*,nu_lw_check,' vs. ',nu_lw(iS)
!                   stop
!                end if
!             end do
!             close(2)

             ! no need to check the pressure layers match
             ! this was done in setup_shortwave
             
!             open(unit=3,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_lw.dat')
!             do iS = 1, nS
!                do iLay = 1, nLay
!                   do iTem = 1,nTem
!                      read(3,'(e12.4)') sigma_lw_ar(iS,iLay,iGas,iTem)
!                   end do
!                end do
!             end do
!             close(3)
             
!          else

             ! calculate cross-sections from scratch
             !   write(6,*) T_k_grid(1,:), T_k_grid(30,:)
             do iLay = 1, size(T_k_grid,1) !nLay
!                if(verbose)then
!                   print*,'For longwave at layer ',iLay,': '
!                endif
                
!                if(nTem==2)then
!                   call line_strength_width(mu_i(iGas),play(iLay),f_i(iLay,iGas)*play(iLay),Tlay(iLay)-dTlay,iGas)
!                   call get_line_abs_cross_section(iGas,nu_lw,sigma_lw,Tlay(iLay)+dTlay)
!                   sigma_lw_ar(:,iLay,iGas,1)    = sigma_lw
                
!                   call line_strength_width(mu_i(iGas),play(iLay),f_i(iLay,iGas)*play(iLay),Tlay(iLay)+dTlay,iGas)
!                   call get_line_abs_cross_section(iGas,nu_lw,sigma_lw,Tlay(iLay)-dTlay)
!                   sigma_lw_ar(:,iLay,iGas,nTem) = sigma_lw
!                else               
                do iTem=1,nTem

                   call line_strength_width(mu_i(iGas),play(1,1,iLay),f_i(1,1,iLay,iGas)*play(1,1,iLay),T_k_grid(iLay,iTem),iGas)
                   call get_line_abs_cross_section(iGas,nu_lw,sigma_lw,T_k_grid(iLay,iTem))
                   sigma_lw_ar(iLay,:,iGas,iTem)    = sigma_lw ! cross-section per molecule of species
!                end if


          ! multiply the species cross-section by molar concentration to get the
          ! contribution to the total cross-section, in cm2/molecule of air.
          ! at this point we also add the CIA and compute the total 'dry' cross-section.
          ! do not add variable gas cross-section to the total yet.
          not_vgas : if( iGas .ne. vgas) then
             if ((f_i(1,1,iLay,iGas)*play(1,1,iLay)>1.0d1)) then
! for CIA calculation
                   sigma_CIA(:) = 0.0d0

                   call get_CIA_cross_section(iGas,nu_lw,sigma_CIA,T_k_grid(iLay,iTem),play(1,1,iLay),f_i(1,1,iLay,:))
                      ! note that f_i is sent as a nGas array to get_CIA_cross_section even when we are inside
                      ! an iGas loop. this is because we need to know the abundance of the partner gas.
                   sigma_total_dry(iLay,:,iTem) = sigma_total_dry(iLay,:,iTem) + &
                      sigma_lw_ar(iLay,:,iGas,iTem)*f_i(1,1,iLay,iGas) + sigma_CIA
             else
                   sigma_total_dry(iLay,:,iTem) = sigma_total_dry(iLay,:,iTem) + &
                      sigma_lw_ar(iLay,:,iGas,iTem)*f_i(1,1,iLay,iGas) 
             end if
          end if not_vgas
                
                   
          end do
          end do

             ! save sigma data for future use
             ! save nu_lw to allow consistency checks
!             open(unit=1,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_lw.dat')
!             open(unit=2,file='saved_sigma_data/nu_lw.dat')
!             do iS = 1, nS
!                do iLay = 1, nLay
!                   do iTem = 1,nTem
!                      write(1,'(e12.4)') sigma_lw_ar(iS,iLay,iGas,iTem)
!                   end do
!                end do
!                write(2,*) nu_lw(iS)
!             end do
!             close(1)
!             close(2)

!          end if
         
       end do big_species_loop
    end if gray_debug_if

  end subroutine setup_longwave

  ! ==================================================================================

  ! ==================================================================================

  subroutine define_cosa_quad(cosa,ang_wt)

    ! calculate the gaussian nodes and weights for
    ! the mean emission angle cosine quadrature

    ! c.f. https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Legendre_quadrature

    implicit none

    real(8), intent(inout) :: cosa(nAng)
    real(8), intent(inout) :: ang_wt(nAng)

    ! gaussian quadrature nodes and weights
    real(8) xi(nAng), ci(nAng)

    ! gauss interval (cosa = 0 to 1 here)
    real(8), parameter :: a = 0.0d0
    real(8), parameter :: b = 1.0d0

    ! move to dynamic allocation eventually

    if(nAng == 1)then
       xi(1) = 0.0d0
       ci(1) = 2.0d0
    elseif(nAng == 2)then
       xi(1) = +sqrt(1.0d0/3.0d0)
       xi(2) = -sqrt(1.0d0/3.0d0)
       ci(:) = 1.0d0
    elseif(nAng == 4)then
       xi(1) = -sqrt(3./7. + (2./7.)*sqrt(6./5.))
       xi(2) = -sqrt(3./7. - (2./7.)*sqrt(6./5.))
       xi(3) = +sqrt(3./7. - (2./7.)*sqrt(6./5.))
       xi(4) = +sqrt(3./7. + (2./7.)*sqrt(6./5.))
       ci(1) = (18.0d0-sqrt(30.0d0))/36.0d0
       ci(2) = (18.0d0+sqrt(30.0d0))/36.0d0
       ci(3) = (18.0d0+sqrt(30.0d0))/36.0d0
       ci(4) = (18.0d0-sqrt(30.0d0))/36.0d0

!!$ elseif(nAng == 8)then
!!$       xi(:) = (/-0.1834346424956498, &
!!$            0.1834346424956498,  &
!!$            -0.5255324099163290, &
!!$            0.5255324099163290,  & 
!!$            -0.7966664774136267, &
!!$            0.7966664774136267,  &
!!$            -0.9602898564975363, &
!!$            0.9602898564975363 /)
!!$       ci(:) =(/ 0.3626837833783620, & 
!!$            0.3626837833783620, &
!!$            0.3137066458778873, &
!!$            0.3137066458778873, &
!!$            0.2223810344533745, &
!!$            0.2223810344533745, &
!!$            0.1012285362903763, &
!!$            0.1012285362903763 /)
    else
       write(*,*) 'Error, emission angle weighting not defined for nAng = ', nAng
    endif

    cosa   = ((b-a)*xi + (b+a))/2
    ang_wt = ci*(b-a)/2
    
    return
  end subroutine define_cosa_quad

  ! ==================================================================================

  subroutine setup_T_k_grid(Tlay,T_k_grid_out)

    real(8), intent(in)  :: Tlay(:)              ! temperature in layer  [K]
    real(8), intent(out) :: T_k_grid_out(:,:) ! temperature array for cross-sections [K]

    integer ierr, iLay, iTem, nLay
    
    nLay = size(Tlay)

    ! read namelist 
!    open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
!    if (ierr/=0) then
!       print*, 'Cannot find required input.nml file, aborting.'
!       call abort
!    else
!       read(10,temperature_k_grid_nml)
!       close(10)
!    endif
    
    ! calculate temperature grid on which cross-sections are defined
    do iLay = 1,nLay
       do iTem = 1,nTem
          T_k_grid_out(iLay,iTem) = Tlay(iLay) + dble(iTem-1-(nTem-1.)/2)*dTlay
       end do
    end do

    ! default values for iT1 and iT2
    iT1(:,:,:) = 1
    iT2(:,:,:) = iT1(:,:,:) + 1

    !T_k_grid_out = T_k_grid

  end subroutine setup_T_k_grid

  ! ==================================================================================

  subroutine calculate_T_k_grid(Tlay,iT1,T_lin_weight)

    real(8), intent(in)  :: Tlay(:,:,:)            ! temperature in layer  [K]
    integer, intent(inout) :: iT1(:,:,:)         ! T-grid array points for linear interpolation [] 
    real(8), intent(out) :: T_lin_weight(:,:,:,:) ! temperature weigting for linear interpolation [] 

    integer :: i, j, iLay, nLay
    
    nLay = size(Tlay,3)

    iT1_out(:,:,:) = -1
    
    ! either just select 1st values in array, or
    ! interpolate over T grid to get cross-section
    ! at actual atmospheric temperature

    ! find lower / upper T values at each layer
do i=1,size(Tlay,1)
   do j=1,size(Tlay,2)

    do iLay = 1,nLay

       find_iT : do 

          if(Tlay(i,j,iLay)<T_k_grid(iLay,iT1(i,j,iLay)))then
             ! move to lower T_k_grid values
             iT1(i,j,iLay) = iT1(i,j,iLay) - 1
             ! halt program if Tlay is below lowest T_k_grid value
             if(iT1(i,j,iLay)==0)then
!                write(*,*) 'Temperature at iLay = ',iLay,' too cold.'
!                write(*,*) 'Tlay(iLay)          = ',Tlay(i,j,iLay),' K.'
!                write(*,*) 'min T_k_grid        = ',T_k_grid(iLay,1),' K.'
!                call error_mesg('radiance: ','too cold for interp',FATAL)
                iT1(i,j,iLay) = 1
                exit
             end if
          elseif(Tlay(i,j,iLay)>T_k_grid(iLay,iT1(i,j,iLay)+1))then
             ! move to higher T_k_grid values
             iT1(i,j,iLay) = iT1(i,j,iLay) + 1
             ! halt program if Tlay is above highest T_k_grid value
             if(iT1(i,j,iLay)==nTem)then
                !write(*,*) 'Temperature at iLay = ',iLay,' too hot.'
                !write(*,*) 'Tlay(iLay)          = ',Tlay(iLay),' K.'
                !write(*,*) 'max T_k_grid        = ',T_k_grid(iLay,nTem),' K.'
!                call error_mesg('radiance: ','too hot for interp',FATAL)
                iT1(i,j,iLay) = nTem-1
                exit
             end if
          else
             ! T_k_grid values are correct, exit
             exit
          end if

       end do find_iT

       T_lin_weight(i,j,iLay,:) = (Tlay(i,j,iLay) - T_k_grid(iLay,iT1(i,j,iLay)))/(T_k_grid(iLay,iT1(i,j,iLay)+1) - T_k_grid(iLay,iT1(i,j,iLay))) 

    end do
   end do
end do
    iT2(:,:,:) = iT1(:,:,:) + 1

!    do iLay = 1,nLay
!       T_lin_weight(:,iLay) = (Tlay(iLay) - T_k_grid(iLay,iT1(iLay)))/(T_k_grid(iLay,iT2(iLay)) - T_k_grid(iLay,iT1(iLay))) 
!    end do
    !if(verbose)then
    !   do iLay = 1,nLay
    !      write(*,'(a15,i3)')    'iT1          = ', iT1(iLay)
    !   end do
    !   do iLay = 1,nLay
    !      write(*,'(a15,f8.3)')  'Tlay         = ', Tlay(iLay)
    !   end do
    !   do iLay = 1,nLay
    !      write(*,'(a15,f8.3)')  'T_lin_weight = ', T_lin_weight(1,iLay)
    !   end do
    !end if

    iT1_out(:,:,:) = iT1

  end subroutine calculate_T_k_grid

  ! ==================================================================================

  subroutine calculate_longwave(iT1,cosa_i,Tlay,Tlev,Ts,T_lin_weight,play,dp,f_i,f_cld,mu_avg, &
       kappa_gray,OLRnu,OLR,I_lev_up_int,I_lev_dn_int)

    !use cross_section,  only : get_CIA_cross_section

    implicit none

    integer, intent(in)  :: iT1(:,:,:)             ! T-grid array points for linear interpolation [] 
    real(8), intent(in)  :: cosa_i                ! emission angle cosine []
    real(8), intent(in)  :: Tlay(:,:,:)            ! temperature in layer  [K]
    real(8), intent(in)  :: Tlev(:,:,:)            ! temperature at level  [K]
    real(8), intent(in)  :: Ts(:,:)                    ! surface temperature [K]
    real(8), intent(in)  :: T_lin_weight(:,:,:,:) ! temperature weigting for linear interpolation [] 
    real(8), intent(in)  :: play(:,:,:)            ! pressure layers [Pa]
    real(8), intent(in)  :: dp(:,:,:)              ! pressure difference [Pa]
    real(8), intent(in)  :: f_i(:,:,:,:)        ! species molar concentration [mol/mol]
    real(8), intent(in)  :: f_cld(:,:,:)        ! cloud mass concentration
    real(8), intent(in)  :: mu_avg(:,:,:)          ! average molar mass in layer [g/mol]
    real(8), intent(in)  :: kappa_gray            ! gray gas mass absorption cross-section for debug [m2/kg]

    real(8), intent(out) :: OLRnu(:,:,:)             ! outgoing longwave spectral irradiance [W/m2/cm^-1/sr] 
    real(8), intent(out) :: OLR(:,:)                   ! total outgoing longwave irradiance [W/m2/sr] 
    real(8), intent(out) :: I_lev_up_int(:,:,:)    ! upward irradiance in each layer [W/m2/sr]
    real(8), intent(out) :: I_lev_dn_int(:,:,:)    ! downward irradiance in each layer [W/m2/sr]

    integer :: i, j, iLay, nLay, iLev, iS, iRev

    nLay = size(Tlay,3)

    sigma_d_i(:,:,:,:) = 0.0d0
    sigma_v_i(:,:,:,:) = 0.0d0
    dtau_lw(:,:,:,:)   = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low
    
    ! ----- beginning of part that requires nTem --------------

    gray_debug_loop : if(gray_debug)then

!       ! nTem = 1 for this case verified in setup_shortwave
!       do i=1,size(Tlay,1)
!       do j=1,size(Tlay,2)
!       do iLay = 1,nLay
!          dtau_lw(i,j,iLay,:) = dtau_lw(i,j,iLay,:)+kappa_gray*dp(i,j,iLay)/grav &
!          * 0.622*f_i(i,j,iLay,1) / (1 - 0.378*f_i(i,j,iLay,1))  !H2O is gray
!       end do
!       end do
!       end do

!       tau_lw_inf(:,:,:) = sum(dtau_lw(:,:,:,:),3)
!       do i=1,size(Tlay,1)
!       do j=1,size(Tlay,2)
!       do iLay = 1,nLay
!          dtau_lw(i,j,iLay,:) = tau_lw_inf(i,j,:)/(1e5)**4 &
!                *4*play(i,j,iLay)**3*dp(i,j,iLay) &
!                + 1.2*dp(i,j,iLay) / 1e5
!       end do
!       end do
!       end do
!       dtau_lw(:,:,:,:) = dtau_lw(:,:,:,:) /2. !ignore cos ang = 0.5 to match
!other fms runs

       do iLay = 1,nLay
          dtau_lw(:,:,iLay,:) = +kappa_gray*dp(1,1,iLay)/grav
       end do

    else

       ! either just select 1st values in array, or
       ! interpolate over T grid to get cross-section
       ! at actual atmospheric temperature
       if(nTem==1)then   
          do i=1,size(Tlay,1)
          do j=1,size(Tlay,2)
             sigma_d_i(i,j,:,:) = sigma_total_dry(:,:,1)
             if(vgas.ne.0) sigma_v_i(i,j,:,:) = sigma_lw_ar(:,:,vgas,1)
          end do
          end do
       else
          iT2(:,:,:) = iT1(:,:,:) + 1

          ! do log-space linear interpolation
          ! calculate logarithm of cross-section arrays
          ! y = y1 + (y2-y1)*(x-x1)/(x2-x1)
          do i=1,size(Tlay,1)
          do j=1,size(Tlay,2)
          do iLay = 1,nLay
             log_sig_d(i,j,iLay,:,:) = log10(sigma_total_dry(iLay,:,iT1(i,j,iLay):iT2(i,j,iLay)) + 1.0d-50)
             if(vgas.ne.0) log_sig_v(i,j,iLay,:,:) = &
             log10(sigma_lw_ar(iLay,:,vgas,iT1(i,j,iLay):iT2(i,j,iLay)) + 1.0d-50)
          end do
          end do
          end do
          sigma_d_i(:,:,:,:) = 10.0d0**(log_sig_d(:,:,:,:,1) + &
          (log_sig_d(:,:,:,:,2) - log_sig_d(:,:,:,:,1))*T_lin_weight(:,:,:,:))
          if(vgas.ne.0) sigma_v_i(:,:,:,:) = 10.0d0**(log_sig_v(:,:,:,:,1) + &
          (log_sig_v(:,:,:,:,2) - log_sig_v(:,:,:,:,1))*T_lin_weight(:,:,:,:))
       end if
       
    
       ! get sigma_CIA for variable gas (continuum absorption if H2O)
       ! and add variable gas cross-section (single-molecule + CIA) to the total
       ! for speed, we assume that _only_ the mixing ratio of the volatile gas
       ! 'iGas_var' may change at each timestep.
       if(vgas==0)then
          sigma_t_i = sigma_d_i
       else
           do i=1,size(Tlay,1)
           do j=1,size(Tlay,2)
           do iLay = 1, nLay

             sigma_CIA(:) = 0.0d0

             ! calculate variable gas continuum at exact instantaneous atmospheric temperature, so no T interpolation
             call get_CIA_cross_section(vgas,nu_lw,sigma_CIA,Tlay(i,j,iLay),play(i,j,iLay),f_i(i,j,iLay,:))
             sigma_t_i(i,j,iLay,:) = sigma_d_i(i,j,iLay,:) + sigma_v_i(i,j,iLay,:)*f_i(i,j,iLay,vgas) + sigma_CIA
          end do
          end do
          end do
       endif

       ! compute total layer optical depth
       ! note array flip [nS,nLay] ==> [nLay,nS]
       do i=1,size(Tlay,1)
       do j=1,size(Tlay,2)
       do iLay = 1, nLay
          dtau_lw(i,j,iLay,:) = dtau_lw(i,j,iLay,:) + &
          1.0d-4*sigma_t_i(i,j,iLay,:)*dp(i,j,iLay)/(grav*mu_avg(i,j,iLay)*mpr)  &
          + 3./4/dens_h2o/r_cld *dp(i,j,iLay)/grav *f_cld(i,j,iLay) !cloud optical thickness assuming absorption cross section of pi*r**2 

       end do
       end do
       end do

    end if gray_debug_loop
    
    ! calculate Planck function vs. height and wavenumber
    do iS = 1, nS
       do i=1,size(Tlay,1)
       do j=1,size(Tlay,2)
       do iLev = 1, nLay+1
          call B_nu(nu_lw(iS),Tlev(i,j,iLev),Bnu_lev(i,j,iLev,iS))
       end do
       do iLay = 1, nLay
          call B_nu(nu_lw(iS),Tlay(i,j,iLay),Bnu_lay(i,j,iLay,iS))
       end do
       call B_nu(nu_lw(iS),Ts(i,j),Bnu_s(i,j,iS))
       end do
       end do
    end do
    
    ! calculate irradiance boundary conditions
    I_lev_up(:,:,1,:)    = Bnu_s(:,:,:)
    I_lev_dn(:,:,nLay+1,:) = 0.0d0

    ! scale dtau_lw by emission angle cosine
    dtau_lw_a(:,:,:,:) = dtau_lw(:,:,:,:)/cosa_i

    ! ----- end of part that requires nTem --------------

    ! vertical path optical depth at TOA
    ! somewhat inefficient to calculate it for every propagation angle
    tau_lw_inf(:,:,:) = sum(dtau_lw(:,:,:,:),3)
    
    ! calculate dTran
    ! the transmission through each layer,
    ! individually, for a given propagation angle.
    dTran = exp(-dtau_lw_a)

    ! make sure transmission does not cause floating point exception
     do i=1,size(Tlay,1)
     do j=1,size(Tlay,2)
     do iLay = 1,nLay
       do iS = 1,nS
          if(dTran(i,j,iLay,iS) < 1.0d-16)then           
             dTran(i,j,iLay,iS) = 1.0d-16
          end if
       end do
    end do
    end do
    end do

    ! calculate Cterm (not dependent on flux direction)
    Cterm(:,:,:,:) = 1.0d0/dtau_lw_a(:,:,:,:) - dTran(:,:,:,:)/(1.0d0 - dTran(:,:,:,:))

    ! Calculate upwards spectral irradiance using Clough et al. (1992) method
    ! starting from the BOA up
    do iLev = 2, nLay+1
       iLay             = iLev - 1
       Beff(:,:,:)      = Bnu_lev(:,:,iLev,:) + &
       2.0d0*(Bnu_lay(:,:,iLay,:) - Bnu_lev(:,:,iLev,:))*Cterm(:,:,iLay,:)
       I_lev_up(:,:,iLev,:) = I_lev_up(:,:,iLev-1,:)*dTran(:,:,iLay,:) + &
       Beff*(1.0d0 - dTran(:,:,iLay,:))
    end do
    ! Calculate downwards spectral irradiance using Clough et al. (1992) method
    ! starting from the TOA down
    do iRev = 1,nLay
       iLev             = nLay+1 - iRev ! start at nLev-1, finish at 1
       iLay             = iLev
       Beff             = Bnu_lev(:,:,iLev,:) + &
       2.0d0*(Bnu_lay(:,:,iLay,:) - Bnu_lev(:,:,iLev,:))*Cterm(:,:,iLay,:)
       I_lev_dn(:,:,iLev,:) = I_lev_dn(:,:,iLev+1,:)*dTran(:,:,iLay,:) + &
       Beff*(1.0d0 - dTran(:,:,iLay,:))
    end do

    ! calculate OLR and irradiances for output
    ! trapezoidal quadrature for the irradiances
    OLRnu(:,:,:)        = I_lev_up(:,:,nLay+1,:) 
    I_lev_up_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
    (I_lev_up(:,:,:,1) + I_lev_up(:,:,:,nS) + 2.0d0*sum(I_lev_up(:,:,:,2:nS-1),4) )
    I_lev_dn_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
    (I_lev_dn(:,:,:,1) + I_lev_dn(:,:,:,nS) + 2.0d0*sum(I_lev_dn(:,:,:,2:nS-1),4) )

!    I_lev_up_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
!    (I_lev_up(:,:,:,2) + I_lev_up(:,:,:,nS-1) + 2.0d0*sum(I_lev_up(:,:,:,2:nS-1),4) )
!    I_lev_dn_int(:,:,:) = ((nu_lw(2)-nu_lw(1))/2.0d0) * &
!    (I_lev_dn(:,:,:,2) + I_lev_dn(:,:,:,nS-1) + 2.0d0*sum(I_lev_dn(:,:,:,2:nS-1),4) )


    OLR(:,:) = I_lev_up_int(:,:,nLay+1)
    !OLR = 0.0d0
    !do iS = 1, nS
    !   OLR = OLR + OLRnu(iS)
    !end do
    !OLR = OLR*dnu_lw

    !-------- end longwave radiative transfer calculation --------

  end subroutine calculate_longwave

  ! ==================================================================================

  subroutine setup_shortwave(T_k_grid,play,ps,f_i,mu_i,nu_sw_out)

    real, intent(in)  :: T_k_grid(:,:) ! temperature array for cross-sections [K] nLay,nTem
    real, intent(in)  :: play(:,:,:)          ! pressure layers [Pa]
    real, intent(in)  :: f_i(:,:,:,:)      ! species molar concentration [mol/mol] nLay,nGas
    real, intent(in)  :: ps(:,:)                  ! surface pressure [Pa]
    real, intent(in)  :: mu_i(nGas)          ! [g/mol]
    real, intent(out) :: nu_sw_out(nS)       ! shortwave wavenumber [cm^-1]    
!    real :: I_temp(nS)                          ! temporary spectral irradiance array [W/m2/cm^-1/sr]
    integer :: nTem_check                          ! number of cross-section temperatures at each layer
    real :: play_check                          ! pressure layer [Pa]
    real :: nu_sw_check                         ! shortwave wavenumber [cm^-1]   
    integer :: iTem, iGas, iLay, iS, i, j, nLay
    nLay = size(play,3)
!    logical, parameter   :: verbose = .true.
    ! read namelist 
!    open(10,file='input.nml',status='old',form='formatted',iostat=ierr)
!    if (ierr/=0) then
!       print*, 'Cannot find required input.nml file, aborting.'
!       call abort
!    else
!       read(10,shortwave_nml)
!       close(10)
!    endif
    ! divide by four to get flux averaged over the sphere
!    ISR = Fstel0/4.0d0
    ! create fine spectral grid and initialize abs array
    nu_sw(1)  = nu_sw1
    dnu_sw = (nu_sw2-nu_sw1)/dble(nS-1)  ! Wavenumber interval [cm^-1]
    do iS = 2, nS
       nu_sw(iS) = nu_sw(iS-1) + dnu_sw
    end do
    nu_sw_out = nu_sw
    ! calculate surface reflectance (albedo)
    Rs(:,:,:) = Asurf
    ! set up the ultraviolet absorbers
!    if(include_UV_abs)then
!       do iGas = 1, nGas
!          if(gas_name(iGas) == 'CO2') call setup_UV('CO2')
!          if(gas_name(iGas) == 'O3_') call setup_UV('O3_')
!       end do
!    end if
    ! note that UV and Rayleigh scattering calculations are currently inconsistent!
    ! calculate stellar Planck spectral irradiance vs. wavenumber
    call init_stellar_spectrum(sigma_sw)
    do i=1,size(play,1)
       do j=1,size(play,2)
          if (cosa0(i,j) .gt. 0) then
             F_stel_dn(i,j,:) = sigma_sw(:)*Fstel0(i,j)/sum(sigma_sw(:)*dnu_sw) ! stellar spectral flux [W/m2/cm^-1]
             I_stel_dn(i,j,:) = F_stel_dn(i,j,:)/(Omega_stel*cosa0(i,j))     ! stellar spectral irradiance [W/m2/cm^-1/sr]
    ! from classic flux defn. (e.g. Goody & Yung)
          else
             F_stel_dn(i,j,:) = 0.
             I_stel_dn(i,j,:) = 0.
          endif
       enddo
    enddo
    ! paint the ground blue with Rayleigh scattering
    ! we can do this when vib-rot bands extend to the near-IR only
    if(.not.multiple_scat .and. .not.rayleigh_top)then 
       call approximate_rayleigh(ps,Rs,Rs_out)
       Rs = Rs_out
       !print*,'HACK: Rayleigh Scattering not added!'
    end if

    if(multiple_scat)then
      !ssa_sw_i(:,:) = 0.0d0 ! initialize single scattering albedo []
      !asy_sw_i(:,:) = 0.0d0 ! initialize asymmetry parameter []
      call setup_mult_scat_core()
    end if 

!    if(verbose)then
!       open(unit=4,file='results/ISRnu.out')
!       do iS = 1, nS
!          write(4,*) F_stel_dn(iS)
!       end do
!       close(4)
!    end if
    !-------- shortwave radiative transfer calculation --------
    ! calculate absorption cross section at each level
    sigma_sw(:)          = 0.0d0
    !tau_sw_inf(:,:,:)    = 0.0d0
    sigma_sw_ar(:,:,:,:) = 0.0d0
    dtau_sw(:,:,:)       = +1.0d-8 ! initialize non-zero to avoid NaN when optical depth is low

    sigma_sw_total_dry(:,:,:) = 0.0d0 ! dry total absorption cross-section (no variable gas) [cm2 / molecules of air] 

    gray_debug_if : if(gray_debug)then
       ! don't calculate anything in this case
!       print *,"gray shortwave"
    else
       big_species_loop : do iGas = 1, nGas
          ! read gas absorption cross section data from file if requested
!          if(.not.calc_sigma(iGas))then
             ! check the number of temperatures for each level match
!             open(unit=1,file='saved_sigma_data/nTem.dat')
!             read(1,*) nTem_check
!             if(nTem_check.ne.nTem)then
!                write(*,*) 'Error: nTem in saved_sigma_data does not match!'
!                stop
!             end if
!             close(1)
             ! check the pressure arrays match
             ! inefficient to do this for each gas but time cost is negligible; whatever
             ! we only do this in shortwave, not in longwave
!             open(unit=2,file='saved_sigma_data/play.dat')
!             do iLay = 1, nLay
!                read(2,*) play_check
!                if(abs(play_check-play(iLay))>1.0d-8)then
!                   write(*,*) 'Error: play data in saved_sigma_data does not match!'
!                   print*,play_check,' vs. ',play(iLay)
!                   stop
!                end if
!             end do
!             close(2)             
             ! check the spectral arrays match
!             open(unit=2,file='saved_sigma_data/nu_sw.dat')
!             do iS = 1, nS
!                read(2,*) nu_sw_check
!                if(abs(nu_sw_check - nu_sw(iS))>1.0d-8)then
!                   write(*,*) 'Error: nu_sw data in saved_sigma_data does not match!'
!                   print*,nu_sw_check,' vs. ',nu_sw(iS)
!                   stop
!                end if
!             end do
!             close(2)             
!             open(unit=3,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_sw.dat')
!             do iS = 1, nS
!                do iLay = 1, nLay
!                   do iTem = 1,nTem
!                      read(3,'(e12.4)') sigma_sw_ar(iS,iLay,iGas,iTem)
!                   end do
!                end do
!             end do
!             close(3)             
!          else
             do iLay = 1,nLay
!                if(verbose)then
!                   print*,'For shortwave at layer ',iLay,': '
!                endif
                ! note no UV-cross-section temperature dependence included for now
!                sigma_UV = 0.0d0
!                if(include_UV_abs)then
!                   if(gas_name(iGas)=='CO2' .or. gas_name(iGas)=='O3_') then
!                      do iS = 1,nS
!                         call calculate_UV(gas_name(iGas),nu_sw(iS),sigma_UV(iS),.false.)
!                      end do
!                   end if
!                end if
                ! calculate cross-sections from scratch
                ! note these line widths will become inaccurate for a gas if
                ! a) it is a major consituent (f_i>0.1 or so) and
                ! b) its abundance varies with time
                do iTem = 1,nTem
                   call line_strength_width(mu_i(iGas),play(1,1,iLay),f_i(1,1,iLay,iGas)*play(1,1,iLay),T_k_grid(iLay,iTem),iGas)
                   call get_line_abs_cross_section(iGas,nu_sw,sigma_sw,T_k_grid(iLay,iTem))
                   sigma_sw_ar(iLay,:,iGas,iTem) = sigma_sw !+ sigma_UV
                   ! cross-section per molecule of species

                   not_vgas : if(.not. iGas == vgas)then
                   if (f_i(1,1,iLay,iGas)*play(1,1,iLay)>1.0d1) then ! include CIA if partial pressure > 10 Pa
                    sigma_CIA(:) = 0.0d0
                    call get_CIA_cross_section(iGas,nu_sw,sigma_CIA,T_k_grid(iLay,iTem),play(1,1,iLay),f_i(1,1,iLay,:))
                    ! note that f_i is sent as a nGas array to get_CIA_cross_section even when we are inside
                    ! an iGas loop. this is because we need to know the abundance of the partner gas.
                    sigma_sw_total_dry(iLay,:,iTem) = sigma_sw_total_dry(iLay,:,iTem) + &
                      sigma_sw_ar(iLay,:,iGas,iTem)*f_i(1,1,iLay,iGas) + sigma_CIA
                   else
                    sigma_sw_total_dry(iLay,:,iTem) = sigma_sw_total_dry(iLay,:,iTem) + &
                      sigma_sw_ar(iLay,:,iGas,iTem)*f_i(1,1,iLay,iGas)

                    end if
                    end if not_vgas
                   
                end do
             end do

             ! save sigma data for future use
             ! save nu_sw, play and nTem to allow consistency checks
!             open(unit=1,file='saved_sigma_data/sigma_'//gas_name(iGas)//'_sw.dat')
!             open(unit=2,file='saved_sigma_data/nu_sw.dat')
!             do iS = 1, nS
!                do iLay = 1, nLay
!                   do iTem = 1,nTem
!                      write(1,'(e12.4)') sigma_sw_ar(iS,iLay,iGas,iTem)
!                   end do
!                end do
!                write(2,*) nu_sw(iS)
!             end do
!             close(1)
!             close(2)
!             open(unit=3,file='saved_sigma_data/nTem.dat')
!             write(3,*) nTem
!             close(3)
!             open(unit=4,file='saved_sigma_data/play.dat')
!             do iLay = 1, nLay
!                write(4,*) play(iLay)
!             end do
!             close(4)
!          end if
       end do big_species_loop
    end if gray_debug_if

  end subroutine setup_shortwave

  ! ==================================================================================
  subroutine init_stellar_spectrum(I_temp)
    real, intent(out) :: I_temp(nS)                ! spectral irradiance array [W/m2/cm^-1/sr]
    integer :: iS_stel, icount, iS
    real :: nu_stel(nS_stel)                          ! stellar spectrum wavenumber [cm^-1]
    real :: F_stel(nS_stel)                           ! stellar spectrum spectral flux (unnormalized) [W/m2/cm^-1]
    I_temp(:) = 0.0d0
    if(stellar_blackbody)then
       do iS = 1, nS
          call B_nu(nu_sw(iS),5800.0d0,I_temp(iS))
       end do
    else
       if(adleo)then
          open(unit=3,file=trim(datadir)//'stellar_data/nu_adleo.dat')
          open(unit=4,file=trim(datadir)//'stellar_data/Fadleo.dat')
       else
          open(unit=3,file=trim(datadir)//'stellar_data/nu.dat')
          open(unit=4,file=trim(datadir)//'stellar_data/Fsol_3p8_Ga.dat')
       end if

       !open(unit=3,file=trim(datadir)//'stellar_data/nu.dat')
       !open(unit=4,file=trim(datadir)//'stellar_data/Fsol_3p8_Ga.dat')
       !open(unit=3,file=trim(datadir)//'stellar_data/nu_adleo.dat')
       !open(unit=4,file=trim(datadir)//'stellar_data/Fadleo.dat')

       ! choice of spectra to be added later (just have filename in input.nml)
       do iS_stel = 1, nS_stel
          read(3,*) nu_stel(iS_stel)
          read(4,*) F_stel(iS_stel)
       end do
       close(3)       
       close(4)
       ! bin the stellar data by band
       ! we assume stellar resolution is higher, and both
       ! wavenumber grids are even
       ! techincally in this code we are using wavenumber bins
       ! cross-sections are evaluated at bin centers
       ! stellar spectrum is added to entire bins
       ! everything is normalized later so units are irrelevant here
       iS     = 1
       icount = 0
       do iS_stel = 1,nS_stel
          if(nu_stel(iS_stel)<nu_sw(iS)-dnu_sw/2.0d0)then
             ! if stellar spectrum point is below starting bin lower boundary, do nothing
          elseif(nu_stel(iS_stel)>nu_sw(iS)+dnu_sw/2.0d0)then
             ! if stellar spectrum point is above current bin upper boundary,
             ! move iS up one bin and add to that, or exit if we were at last bin.
             if(iS<nS)then
                I_temp(iS) = I_temp(iS)/dble(icount)      ! normalise amount in previous bin
                iS         = iS + 1                       ! move up one bin
                I_temp(iS) = I_temp(iS) + F_stel(iS_stel) ! add to total in new bin
                icount     = 1                            ! icount is reset
             else
                exit
             end if
          else
             ! if stellar spectrum point falls in current bin, just add to total.
             I_temp(iS) = I_temp(iS) + F_stel(iS_stel)
             icount     = icount + 1 ! one point added to current bin
          end if
       end do
    end if
    return
  end subroutine init_stellar_spectrum
  
  ! ==================================================================================

  subroutine approximate_rayleigh(ps,Rbelow,Rabove)
    real, intent(in)  :: ps(:,:)         ! surface pressure [Pa]
    real, intent(in)  :: Rbelow(:,:,: ) ! reflectivity/albedo below Rayleigh layer []
    real, intent(out) :: Rabove(:,:,: ) ! reflectivity/albedo above Rayleigh layer []
    ! no temperature dependence so we only need to do this once
    ! unless we have really variable major constituents in the atmosphere

    ! general variables
    real(8), parameter :: gamma = 1d0 !0.5d0*sqrt(3.0d0) ! hemi-isotropic approximation
    !real, parameter   :: gamma = 3.0d0/4.0d0       ! Eddington approximation
    real, dimension(size(Rbelow,1), size(Rbelow,2), size(Rbelow,3)) :: beta, f, R_dn, R_up, tau_Ray
    real, dimension(size(Rbelow,3)) :: lam_sw
    integer :: i, j
!    real(8) beta(nS)
!    real(8) f(nS)
!    real(8) R_dn(nS)    ! reflectivity for downward beam
!    real(8) R_up(nS)    ! reflectivity for upward beam
    !real(8) Rpla(nS)    ! planetary (net) reflectivity
!    real(8) tau_Ray(nS) ! total vertical path Rayleigh optical depth []
!    real(8) lam_sw(nS)  ! wavelength [um]
    
    ! asymmetry parameter is zero for Rayleigh scattering
    lam_sw  = 1.0d4/nu_sw 
    ! kg/m2 x m2/kg = []
    ! At 1.5 bar CO2 on Mars albedo is 0.02 greater for 2nd formula
    !tau_Ray = (ps/grav)*1.2d-6*lam_sw**(-4) ! PPC Table 5.2
    do i=1,size(ps,1)
       do j=1,size(ps,2)
          if (cosa0(i,j) .gt. 0) then
             tau_Ray(i,j,:) = sigma_sca *lam_sw(:)**(-4)*(ps(i,j)/grav)
!             tau_Ray(i,j,:) = sigma_sca *lam_sw(:)**(-4)*(1.0d0 + 0.013d0*lam_sw(:)**(-2))*(ps(i,j)/grav)          ! Hansen & Travis eqn. (2.32)
             beta(i,j,:)    = 1.0d0 - exp(-tau_Ray(i,j,:)/cosa0(i,j))
             f(i,j,:)       = gamma*tau_Ray(i,j,:)
             R_dn(i,j,:)    = ((0.5d0 - gamma*cosa0(i,j) )*beta(i,j,:) + f(i,j,:) )/(1.0d0 + f(i,j,:) )                             ! downwards beam albedo
             R_up(i,j,:)    = f(i,j,:) /(1.0d0 + f(i,j,:) )                                                            ! upwards beam albedo
          else
             R_dn(i,j,:) = 0.
             R_up(i,j,:) = 0.
          endif
       enddo
    enddo
    !Rpla  = 1.0d0 - (1 - R_up)*(1 - Rs)/((1.0d0 - Rs)*R_dn + (1.0d0 - R_dn)) ! planetary albedo
    Rabove  = 1.0d0 - (1 - R_dn)*(1 - Rbelow)/((1.0d0 - Rbelow)*R_up + (1.0d0 - R_up)) ! planetary albedo
    ! should work for clear-sky atmospheres because Rayleigh scattering and
    ! molecular absorption occur in separate regions.
    ! literally we are painting the ground (or TOA) blue.
    !if(rayleigh_top)then ! Rayleigh scattering applied at TOA
    !   
    !else                 ! Rayleigh scattering applied at surface 
    !   Rs = Rpla
    !endif
  end subroutine approximate_rayleigh

  ! ==================================================================================

  subroutine calculate_shortwave(iT1,cosa_i,T_lin_weight,ts,dp,f_i,f_cld,mu_avg,kappa_gray,bucket,&
       I_up_TOA_nu,I_up_TOA,I_dn_BOA,I_lev_up_int,I_lev_dn_int,I_lev_dir_int)
    
    integer, intent(in)  :: iT1(:,:,:)             ! T-grid array points for linear interpolation []nLay 
    real(8), intent(in)  :: cosa_i                ! emission angle cosine []
    real(8), intent(in)  :: T_lin_weight(:,:,:,:) ! temperature weigting for linear interpolation []nS,nLay 
    real(8), intent(in)  :: ts(:,:)                    ! surface temperature [K]
    real(8), intent(in)  :: dp(:,:,:)              ! pressure difference [Pa]
    real(8), intent(in)  :: f_i(:,:,:,:)        ! species molar concentration [mol/mol]
    real(8), intent(in)  :: f_cld(:,:,:)        ! cloud mass concentration
    real(8), intent(in)  :: mu_avg(:,:,:)          ! average molar mass in layer [g/mol]
    real(8), intent(in)  :: kappa_gray            ! gray gas mass absorption cross-section for debug [m2/kg]
    real(8), intent(in)  :: bucket(:,:)           !bucket water depth [m]
    
    real(8), intent(out) :: I_dn_BOA(:,:)               ! ground diffuse irradiance [W/m2/sr] 
    real(8), intent(out) :: I_up_TOA_nu(:,:,:)        ! outgoing diffuse spectral irradiance [W/m2/cm^-1/sr] 
    real(8), intent(out) :: I_up_TOA(:,:)               ! outgoing diffuse irradiance [W/m2/sr] 
    real(8), intent(out) :: I_lev_up_int(:,:,:)     ! upward irradiance in each layer [W/m2/sr]
    real(8), intent(out) :: I_lev_dn_int(:,:,:)     ! downward irradiance at each level [W/m2/sr]
    real(8), intent(out) :: I_lev_dir_int(:,:,:)    ! downward direct irradiance at each level [W/m2/sr], now converted to direct beam flux

    real(8),dimension(size(I_up_TOA_nu,1),size(I_up_TOA_nu,2),size(I_up_TOA_nu,3)) ::  I_dn_BOA_nu
    ! ground diffuse spectral irradiance [W/m2/cm^-1/sr]
    !real, dimension(size(I_lev_up_int,1), size(I_lev_up_int,2), size(I_lev_up_int,3)) :: I_lev_dn_int                    ! downward (direct only) irradiance in each layer [W/m2/sr]
    integer :: i,j, iLay, nLay, iLev, nLev, iRev
    integer :: iGas,iS
    nLay = size(dp, 3)
    nLev = size(I_lev_up_int, 3)
    ! note: we are currently neglecting H2O continuum effects in the shortwave
    ! fine for cold climates, will cause problems for RG calculations etc.
    ! as a result dtau_sw calculation is done a little differently from that for dtau_lw
    ! it is not at all optimized for speed; this will be tackled later.
    ! ----- beginning of part that requires nTem --------------
    dtau_sw_i(:,:,:,:)   = 1.0d-18
    dtau_sw_ray(:,:,:,:) = 1.0d-18
    dtau_sw_cld(:,:,:,:) = 1.0d-18
    gray_debug_loop: if(gray_debug)then
       if(nTem>1)then
          write(*,*) 'nTem>1 and gray_debug incompatible.'
          stop
       end if
!no shortwave abs
       do iS = 1,nS
          dtau_sw_i(:,:,:,iS) = +kappa_gray*dp(:,:,:)/grav
       end do

    else

       ! either just select 1st values in array, or
       ! interpolate over T grid to get cross-section
       ! at actual atmospheric temperature
       if(nTem==1)then   
          do i=1,size(Tlay,1)
          do j=1,size(Tlay,2)
             sigma_d_i(i,j,:,:) = sigma_sw_total_dry(:,:,1)
             if(vgas.ne.0) sigma_v_i(i,j,:,:) = sigma_sw_ar(:,:,vgas,1)
          end do
          end do
       else
          iT2(:,:,:) = iT1(:,:,:) + 1

          ! do log-space linear interpolation
          ! calculate logarithm of cross-section arrays
          ! y = y1 + (y2-y1)*(x-x1)/(x2-x1)
          do i=1,size(Tlay,1)
          do j=1,size(Tlay,2)
          do iLay = 1,nLay
             log_sig_d(i,j,iLay,:,:) = log10(sigma_sw_total_dry(iLay,:,iT1(i,j,iLay):iT2(i,j,iLay)) + 1.0d-50)
             if(vgas.ne.0) log_sig_v(i,j,iLay,:,:) = &
             log10(sigma_sw_ar(iLay,:,vgas,iT1(i,j,iLay):iT2(i,j,iLay)) + 1.0d-50)
          end do
          end do
          end do
          sigma_d_i(:,:,:,:) = 10.0d0**(log_sig_d(:,:,:,:,1) + &
          (log_sig_d(:,:,:,:,2) - log_sig_d(:,:,:,:,1))*T_lin_weight(:,:,:,:))
          if(vgas.ne.0) sigma_v_i(:,:,:,:) = 10.0d0**(log_sig_v(:,:,:,:,1) + &
          (log_sig_v(:,:,:,:,2) - log_sig_v(:,:,:,:,1))*T_lin_weight(:,:,:,:))
       end if
       
    
       ! get sigma_CIA for variable gas (continuum absorption if H2O)
       ! and add variable gas cross-section (single-molecule + CIA) to the total
       ! for speed, we assume that _only_ the mixing ratio of the volatile gas
       ! 'iGas_var' may change at each timestep.
       if(vgas==0)then
          sigma_t_i = sigma_d_i
       else
           do i=1,size(Tlay,1)
           do j=1,size(Tlay,2)
           do iLay = 1, nLay

             sigma_CIA(:) = 0.0d0

             ! calculate variable gas continuum at exact instantaneous
             ! atmospheric temperature, so no T interpolation
             call get_CIA_cross_section(vgas,nu_sw,sigma_CIA,Tlay(i,j,iLay),play(i,j,iLay),f_i(i,j,iLay,:))
             sigma_t_i(i,j,iLay,:) = sigma_d_i(i,j,iLay,:) +sigma_v_i(i,j,iLay,:)*f_i(i,j,iLay,vgas) + sigma_CIA
          end do
          end do
          end do
       endif

       ! compute total layer optical depth
       ! note array flip [nS,nLay] ==> [nLay,nS]
       do i=1,size(Tlay,1)
       do j=1,size(Tlay,2)
       do iLay = 1, nLay
          dtau_sw_i(i,j,iLay,:) = dtau_sw_i(i,j,iLay,:) + &
          1.0d-4*sigma_t_i(i,j,iLay,:)*dp(i,j,iLay)/(grav*mu_avg(i,j,iLay)*mpr)
          !rayleigh scattering optical thickness of Earth's air at 1um 
          dtau_sw_ray(i,j,iLay,:) = sigma_sca *dp(i,j,iLay)/grav *(1e4/nu_sw(:))**(-4)
          !      0.3066*2.49e-6 *dp(i,j,iLay)/grav *&
          !      (1e4/nu_sw(:))**(-4)
          !cloud scattering optical thickness 
          dtau_sw_cld(i,j,iLay,:) = 1.5/dens_h2o/r_cld *dp(i,j,iLay)/grav *f_cld(i,j,iLay)
          !dtau_sw_cld(i,j,iLay,:) = 1e-10 
       end do
       end do
       end do

    end if gray_debug_loop

    ! ----- end of part that requires nTem --------------
    ! dtau_sw_i is now the optical thickness by atmospheric absorption
    dtau_sw_i = dtau_sw_i + dtau_sw_ray + dtau_sw_cld !now the total

    ! calculate vertical path total optical depth at TOA
    ! somewhat inefficient to calculate it for every propagation angle
    !tau_sw_inf(:,:,:) = sum(dtau_sw_i(:,:,:,:),3) 
    ! calculate irradiance boundary conditions
    !I_lev_dir(:,:,nLev,:) = I_stel_dn(:,:,:)
    !do i=1,size(dp,1)
    !do j=1,size(dp,2)
    !   if (cosa0(i,j) .gt. 0) then
    ! calculate transmission
    ! note this is vertical path optical depth here
    !      dTran0(i,j,:,:) = exp(-dtau_sw_i(i,j,:,:)/cosa0(i,j))  ! layer transmission for incoming stellar radiation (direct beam)
    !   else
    !      dTran0(i,j,:,:) = 0.
    !   endif
    !enddo
    !enddo
    ! Calculate downwards (direct beam) spectral irradiance 
    ! starting from the TOA down
    !do iRev = 1,nLay
    !   iLev                  = nLev - iRev ! start at nLev-1, finish at 1
    !   iLay                  = iLev
    !   I_lev_dir(:,:,iLev,:) = I_lev_dir(:,:,iLev+1,:)*dTran0(:,:,iLay,:) 
    !end do
    ! --- up to this point, nothing depends on emission angle ---
    ! --- make this into a separate subroutine eventually ---
    !dTran(:,:,:,:)  = exp(-dtau_sw_i(:,:,:,:)/cosa_i) ! layer transmission for radiation scattered from surface
    ! make sure transmission does not cause floating point exception
!    do iLay = 1,nLay
!       do iS = 1,nS
!          if(dTran(i,j,iLay,iS) < 1.0d-16)then           
!             dTran(i,j,iLay,iS) = 1.0d-16
!          end if
!       end do
!    end do
    where (dTran .lt. 1.0d-16) dTran = 1.0d-16

    ! now calculate the diffuse upward and downward irradiances
    if(multiple_scat)then
        ! eventually, use these:
        ! ssa_sw_i(nLay,nS) ! single scattering albedo of each layer at local
        ! temperature []
        ! asy_sw_i(nLay,nS) ! asymmetry parameter of each layer at local
        ! temperature []
        !aL(:,:,:) = 0.9999 !0.4 !d-16 ! single scattering albedo []
        !bL(:,:,:) = bl_cld !0.5d0 !0.5d0*(1-gsca*3*0.5**2) ! backscatter parameter []
        ! g = -1, b = 1.0: complete backscattering
        ! g =  0, b = 0.5: isotropic scattering
        ! g = +1, b = 0.0: complete forward scattering
       do iS = 1,nS
        ! CAREFUL WITH I_stel_dn vs. F0
        ! set muav somewhere
        ! create correct aL, bL and tauL arrays 
        ! flip tau to go from PPC convention to Thomas et al. (standard)
        ! convention
        ! with tau=0 at TOA and tau=tau_max at BOA.
        ! it'd be nice to avoid this, but such is life.

    ssa_sw_i  = ( dtau_sw_ray(:,:,:,iS) + dtau_sw_cld(:,:,:,iS)) /&
            dtau_sw_i(:,:,:,iS)
    delta_f     = dtau_sw_cld(:,:,:,iS)/&
            ( dtau_sw_ray(:,:,:,iS) + dtau_sw_cld(:,:,:,iS)) * gsca**2
    x1_rescaled = dtau_sw_cld(:,:,:,iS)/&
            ( dtau_sw_ray(:,:,:,iS) + dtau_sw_cld(:,:,:,iS)) *gsca !before rescaling
    x1_rescaled = (x1_rescaled - delta_f) / (1.-delta_f)  !delta-two-term scaling
    !aL          = (1.-delta_f) /(1.-ssa_sw_i*delta_f) *ssa_sw_i
    !bL          = 1./2 - 3./8*x1_rescaled
    
!aL, bL should also be flipped
        do iLay=1,nLay
           iRev=nLay-iLay+1
           aL(:,:,iLay) = (1.-delta_f(:,:,iRev)) /&
                   (1.-ssa_sw_i(:,:,iRev)*delta_f(:,:,iRev)) *&
                   ssa_sw_i(:,:,iRev)
           bL(:,:,iLay) = 1./2 - 3./8*x1_rescaled(:,:,iRev) 
           b0L(:,:,iLay)= 1./2*(1.-3./2*x1_rescaled(:,:,iRev)*cosa0(:,:)) !back scattering parameter for direct beam
        enddo        

        where (aL .gt. 1.-1e-5) aL=1.-1e-5
        where (aL .lt. 1e-5)    aL=1e-5
        !aL = 1e-5
        !if (iS == 188) print *, dtau_sw_i(1,1,:,iS)
        !print *, aL(1,1,:)

        tauL(:,:,1) = 0.0d0
        do iLev=2,nLev
           iRev = nLev-iLev+1
           tauL(:,:,iLev) = tauL(:,:,iLev-1) + dtau_sw_i(:,:,iRev,iS) *&
                   (1.-ssa_sw_i(:,:,iRev)*delta_f(:,:,iRev))
!           tauL(iLev) = dble(iLev-1)
        end do

        ! calculate the direct and diffuse fluxes
        do i=1,size(dp,1)
           do j=1,size(dp,2)
              if (cosa0(i,j) .gt. 0) then !dayside hemisphere
                 I_lev_dir_ms(i,j,:) = I_stel_dn(i,j,iS) &
                         * exp(-tauL(i,j,:) /cosa0(i,j))

                 if (bucket(i,j) .le. 0) then
                    call get_Iup_Idn(aL(i,j,:),bL(i,j,:),b0L(i,j,:),tauL(i,j,:),&
                I_stel_dn(i,j,iS)*(Omega_stel),Asurf,cosa_i,cosa0(i,j),&
                I_lev_up_ms(i,j,:),I_lev_dn_ms(i,j,:),.false.,nLay,nLev)
 
                 else if (ts(i,j) .gt. tfreeze) then 
                 
                    call get_Iup_Idn(aL(i,j,:),bL(i,j,:),b0L(i,j,:),tauL(i,j,:),&
                I_stel_dn(i,j,iS)*(Omega_stel),A_water,cosa_i,cosa0(i,j),&
                I_lev_up_ms(i,j,:),I_lev_dn_ms(i,j,:),.false.,nLay,nLev)

                 else
                    call get_Iup_Idn(aL(i,j,:),bL(i,j,:),b0L(i,j,:),tauL(i,j,:),&
                I_stel_dn(i,j,iS)*(Omega_stel),A_ice,cosa_i,cosa0(i,j),&
                I_lev_up_ms(i,j,:),I_lev_dn_ms(i,j,:),.false.,nLay,nLev)
                 endif

              else
                 I_lev_dir_ms(i,j,:)= 0.
                 I_lev_up_ms(i,j,:) = 0.
                 I_lev_dn_ms(i,j,:) = 0.
              endif
           enddo
        enddo
        ! flip the fluxes to go from Thomas et al. to PPC convention
        do iLev=1,nLev
           iRev = nLev-iLev+1
           I_lev_up(:,:,iLev,iS) = I_lev_up_ms(:,:,iRev)
           I_lev_dn(:,:,iLev,iS) = I_lev_dn_ms(:,:,iRev)
           I_lev_dir(:,:,iLev,iS)= I_lev_dir_ms(:,:,iRev)
        end do

       end do 

!       print *, I_lev_up_ms(1,1,:), I_lev_dn_ms(1,1,:)

     else 

      ! Assume downward diffuse spectral irradiance at TOA is zero
      I_lev_dn(:,:,:,:) = 0.0d0 

    ! Calculate upwards spectral irradiance 
    ! starting from the BOA up
    ! assuming a Lambertian surface
    ! F_dn = Omega_stel*I_dn*cosa0
    ! upwards is isotropic so by defn. I not a fn. of angle
    ! then F_up = 2*pi*I_up*sum(cosa_i*ang_wt) = pi*I_up
    ! conservation of energy: F_up = F_dn*Rs
    ! so  pi*I_up = Omega_stel*I_dn*cosa0*Rs
    do i=1,size(dp,1)
    do j=1,size(dp,2)
       I_lev_up(i,j,1,:) = I_lev_dn(i,j,1,:)*Rs(i,j,:)*(Omega_stel*cosa0(i,j)/pi)
    enddo
    enddo

    do iLev = 2, nLev
       iLay             = iLev - 1
       I_lev_up(:,:,iLev,:) = I_lev_up(:,:,iLev-1,:)*dTran(:,:,iLay,:)
    end do
    endif !end multiple scattering 

!    I_lev_up_int(:,:,:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_lev_up(:,:,:,1) + I_lev_up(:,:,:,nS) + 2.0d0*sum(I_lev_up(:,:,:,2:nS-1),4) )
!    I_lev_dn_int(:,:,:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_lev_dn(:,:,:,1) + I_lev_dn(:,:,:,nS) + 2.0d0*sum(I_lev_dn(:,:,:,2:nS-1),4) )
!    GSR(:,:) = F_lev_dn_dir(:,:,1)
!    OSR(:,:) = I_lev_up_int(:,:,nLev)

    ! GSR and OSR are a bit of a mess now
    ! problem is that you can no longer separate direct + diffuse
    ! and you want to preserve ability to do multi-stream avg in some cases

    ! solution: add direct part after

    ! calculate spectral fluxes and averaged irradiances for output
    ! trapezoidal quadrature for the irradiances
    ! CHECK THIS 
    ! todo: replace pi, Omega_stel*cosa0 with S_dif, S_dir
    !OSRnu(:)        = I_lev_up(nLev,:)               ! [W/m2/cm^-1/sr]
    !(notation not ideal! this is an irradiance)
    !GSRnu(:)        = I_lev_dn(1,:) !+ Omega_stel*cosa0*I_lev_dir(1,:) !
    ![W/m2/cm^-1/sr]
    I_dn_BOA_nu(:,:,:)        = I_lev_dn(:,:,1,:)    ! diffuse radiation at bottom of atmosphere (used to get GSR) [W/m2/cm^-1/sr], 1 is surface now
    I_up_TOA_nu(:,:,:)        = I_lev_up(:,:,nLev,:) ! diffuse radiation at top of atmosphere (used to get OSR) [W/m2/cm^-1/sr]

    ! if required, calculate Rayleigh scattering at TOA
    if((.not. multiple_scat) .and. rayleigh_top)then
      call error_mesg('radiance: ','This option is no longer permitted, or needed',FATAL)
    end if

    ! sum of diffuse and direct downward beams 
    ! re-check differences between trapezoidal and riemann integration
    !I_lev_up_int(:)  = dnu_sw*sum(I_lev_up(:,1:nS), 2) 
    !I_lev_dn_int(:)  = dnu_sw*sum(I_lev_dn(:,1:nS), 2) 
    !I_lev_dir_int(:) = dnu_sw*sum(I_lev_dir(:,1:nS),2) 

    I_lev_up_int(:,:,:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_lev_up(:,:,:,1) + &
        I_lev_up(:,:,:,nS) + 2.0d0*sum(I_lev_up(:,:,:,2:nS-1),4) )
    I_lev_dn_int(:,:,:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_lev_dn(:,:,:,1) + &
        I_lev_dn(:,:,:,nS) + 2.0d0*sum(I_lev_dn(:,:,:,2:nS-1),4) )
    I_lev_dir_int(:,:,:)= ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_lev_dir(:,:,:,1) + &
        I_lev_dir(:,:,:,nS) + 2.0d0*sum(I_lev_dir(:,:,:,2:nS-1),4) )

    I_dn_BOA(:,:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_dn_BOA_nu(:,:,1) + &
        I_dn_BOA_nu(:,:,nS) + 2.0d0*sum(I_dn_BOA_nu(:,:,2:nS-1),3) )
    I_up_TOA(:,:) = ((nu_sw(2)-nu_sw(1))/2.0d0) * (I_up_TOA_nu(:,:,1) + &
        I_up_TOA_nu(:,:,nS) + 2.0d0*sum(I_up_TOA_nu(:,:,2:nS-1),3) )
    
    do iLev=1,nLev
       I_lev_dir_int(:,:,iLev) =  I_lev_dir_int(:,:,iLev)*omega_stel*cosa0(:,:) 
    enddo
    ! note we could have done whole direct beam calculation with fluxes
    ! (Omega_stel value is irrelevant).
    !I_dn_BOA(:,:) = 0.0d0
    !I_up_TOA(:,:) = 0.0d0
    !do iS = 1, nS
    !   I_dn_BOA = I_dn_BOA + I_dn_BOA_nu(iS)
    !   I_up_TOA = I_up_TOA + I_up_TOA_nu(iS)
    !end do
    !I_dn_BOA = I_dn_BOA*dnu_sw
    !I_up_TOA = I_up_TOA*dnu_sw

!---------TEMPORARY----------------------------
!    open(unit=5,file='results/Ilev_dn_sw_int.out')
!    open(unit=6,file='results/Ilev_up_sw_int.out')
!    open(unit=7,file='results/Ilev_dir_sw_int.out')
!       do iLev = 1, nLev
!          if(I_lev_dn_int(iLev)<1.0d-30) I_lev_dn_int(iLev) = 0.0d0
!          if(I_lev_up_int(iLev)<1.0d-30) I_lev_up_int(iLev) = 0.0d0
!          if(I_lev_dir_int(iLev)<1.0d-30) I_lev_dir_int(iLev) = 0.0d0
!          ! set very small values to zero so that saved quantities are not
!          ! messed up
!          write(5,'(e12.4)') I_lev_dn_int(iLev)
!          write(6,'(e12.4)') I_lev_up_int(iLev)
!          write(7,'(e12.4)') I_lev_dir_int(iLev)
!       end do
!    close(5)
!    close(6)
!    close(7)
! --------------------------------------------

    !-------- end shortwave radiative transfer calculation -----


  end subroutine calculate_shortwave

  ! ==================================================================================

subroutine orbital_dis(kappa_season, distance)
!real(8), intent(in), dimension(:) :: t, q, p_full
real(8), intent(in) :: kappa_season
real(8), intent(out):: distance

distance = semimajor*(1.-ecc*ecc)/(1.+ecc*cos(kappa_season))
end subroutine orbital_dis

subroutine kappadot(kappa_season, vel)
real(8), intent(in) :: kappa_season
real(8), intent(out):: vel

real(8) :: semiminor, dis
semiminor = semimajor *(1-ecc**2)**0.5
call orbital_dis(kappa_season, dis)
vel = semimajor * semiminor * 2* pi /Torb &
        /seconds_per_day /dis**2
end subroutine kappadot

subroutine find_substellar(kappa_old, lon_old, &
        kappa_new, lon_new, lat_new)
real(8), intent(in) :: kappa_old, lon_old
real(8), intent(out):: kappa_new, lon_new, lat_new

real(8) :: kappavel
call kappadot(kappa_old, kappavel)

kappa_new = kappa_old + kappavel* deltat_rad
lat_new   = asin( sin(obl*deg_to_rad)*cos(kappa_new -prec*deg_to_rad))
lon_new   = lon_old - 2*pi*(1./Trot - 1./Torb) &
        /seconds_per_day *deltat_rad

if (lon_new < 0) lon_new = lon_new + 2*pi

end subroutine find_substellar


end module radiance_mod





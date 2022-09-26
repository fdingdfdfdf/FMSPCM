!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for FV dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------


use constants_mod, only: grav, kappa, cp_air, pi, rdgas, rvgas, SECONDS_PER_DAY, &
                         cp_vapor, cp_water, latent_heat, rho_cp, rho0, r_vir, cp_vir, &
                         ps0, deltat_rad, dens_h2o, hlv

use fms_mod,       only: file_exist, open_namelist_file,   &
                         error_mesg, FATAL,                &
                         check_nml_error, stdlog, stdout,  &
                         write_version_number,             &
                         close_file, set_domain, nullify_domain, mpp_pe, mpp_root_pe, &
                         read_data, field_exist
use time_manager_mod, only: time_type, get_time, set_time, operator(+)
use mpp_domains_mod,  only: domain2d

!------------------
! FV specific codes:
!------------------
use fv_arrays_mod, only: fv_atmos_type
use fv_control_mod,only: fv_init, domain, fv_end, adiabatic, p_ref
use fv_phys_mod,   only: fv_phys, fv_nudge
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time
use fv_timing_mod,   only: timing_on, timing_off
use fv_restart_mod, only: fv_restart
use fv_dynamics_mod, only: fv_dynamics
use fv_grid_tools_mod, only: grid_type, globalsum, area
!use lin_cld_microphys_mod, only: lin_cld_microphys_init, lin_cld_microphys_end

!FD, LBL
use diag_manager_mod,   only: register_diag_field, send_data
use radiance_mod,       only: radiance_init, radiance_update
!use qe_moist_convection_mod, only: dry_convection 
!FD, PBL
use vert_turb_driver_mod, only: vert_turb_driver_init, vert_turb_driver, vert_turb_driver_end
use vert_diff_mod,      only: vert_diff_init, gcm_vert_diff_down, gcm_vert_diff_up, vert_diff_end, surf_diff_type
use surface_flux_mod,   only: surface_flux
use mixed_layer_mod,    only: mixed_layer_init, mixed_layer, mixed_layer_end
use fv_update_phys_mod, only: update_dwinds_phys 
use fv_mp_mod,          only: gid, masterproc
!FD, moist
use lscale_cond_mod, only: lscale_cond_init, lscale_cond
use  qe_moist_convection_mod, only: qe_moist_convection_init, &
                                    qe_moist_convection,   &
                                    qe_moist_convection_end
use manabe_convection_mod, only: moist_convection
use fv_grid_utils_mod,     only: g_sum
use fms_io_mod,            only: get_tile_string, field_size
use cloud_evap_mod,        only: calculate_rehum, cloud_evap
!-----------------------------------------------------------------------

implicit none
private

public   atmosphere_init, atmosphere,  atmosphere_end, atmosphere_domain

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 18.0 2010/03/02 23:26:58 fms Exp $'
character(len=128) :: tag = '$Name: riga $'
character(len=10), parameter :: mod_name='atmosphere'

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos
real                    :: dt_atmos
integer :: sec
integer days, seconds

logical :: cold_start      = .false.       ! read in initial condition

type(fv_atmos_type), allocatable :: Atm(:)
!-----------------------------------------------------------------------
!FD
logical :: do_manabe      = .true.
logical :: do_bettsmiller = .false.
!FD, PBL
logical :: used, doing_edt, doing_entrain
logical :: module_is_initialized =.false.
logical :: turb = .false.
logical :: do_virtual = .false. ! whether virtual temp used in gcm_vert_diff
real :: roughness_heat = 0.05
real :: roughness_moist = 0.05
real :: roughness_mom = 0.05
! bucket model,adpated from Schneider's FMS
logical :: bucket = .true. 
real :: init_bucket_depth = 0.01 ! default initial bucket depth in m LJJ
real :: init_bucket_depth_land = 20. 
real :: land_left       = 90.
real :: land_right      = 270.
real :: land_bottom     = -10.
real :: land_top        = 10.
real :: max_bucket_depth_land = 1000. ! default large value
real :: robert_bucket = 0.04   ! default robert coefficient for bucket depth LJJ
real :: raw_bucket = 0.53       ! default raw coefficient for bucket depth LJJ
real :: damping_coeff_bucket = 200.
! default damping coefficient for diffusing of surface water - default is no diffusion. [degrees/year] 
! end addition by LJJ
real :: bucket_scale = 100 !in yrs
real :: bucket_i, bucket_f
logical :: read_bucket = .false.
logical :: bucket_iter = .true.
logical :: bucket_iter_res = .true.

real :: cloud_time = 3600.

namelist/atmosphere_nml/ do_bettsmiller, do_manabe, bucket, init_bucket_depth, &
        bucket_iter, bucket_iter_res, bucket_scale, &
        cloud_time

!===================================================================
real, allocatable, dimension(:,:)   ::                                        &
     z_surf,               &   ! surface height
     t_surf,               &   ! surface temperature
     q_surf,               &   ! surface moisture
     u_surf,               &   ! surface U wind
     v_surf,               &   ! surface V wind
     rough_mom,            &   ! momentum roughness length for surface_flux
     rough_heat,           &   ! heat roughness length for surface_flux
     rough_moist,          &   ! moisture roughness length for surface_flux
     gust,                 &   ! gustiness constant
     z_pbl,                &   ! gustiness constant
     flux_t,               &   ! surface sensible heat flux
     flux_q,               &   ! surface moisture flux
     flux_r,               &   ! surface radiation flux
     flux_u,               &   ! surface flux of zonal mom.
     flux_v,               &   ! surface flux of meridional mom.
     drag_m,               &   ! momentum drag coefficient
     drag_t,               &   ! heat drag coefficient
     drag_q,               &   ! moisture drag coefficient
     w_atm,                &   ! wind speed
     ustar,                &   ! friction velocity
     bstar,                &   ! buoyancy scale
     qstar,                &   ! moisture scale
     dhdt_surf,            &   ! d(sensible heat flux)/d(surface temp)
     dedt_surf,            &   ! d(latent heat flux)/d(surface temp)???
     dedq_surf,            &   ! d(latent heat flux)/d(surface moisture)???
     drdt_surf,            &   ! d(upward longwave)/d(surface temp)
     dhdt_atm,             &   ! d(sensible heat flux)/d(atmos.temp)
     dedq_atm,             &   ! d(latent heat flux)/d(atmospheric mixing rat.)
     dtaudv_atm,           &   ! d(stress component)/d(atmos wind)
     dtaudu_atm,           &   ! d(stress component)/d(atmos wind)
     fracland,             &   ! fraction of land in gridbox
     rough                     ! roughness for vert_turb_driver

real, allocatable, dimension(:,:,:) :: tg_tmp, qg_tmp, conv_cloud, rain_profile
real, allocatable, dimension(:,:,:) :: p_full, p_half, z_full, z_half
real, allocatable, dimension(:,:,:) :: dt_ug, dt_vg, dt_tg
real, allocatable, dimension(:,:,:) :: dt_ua, dt_va
real, allocatable, dimension(:,:,:,:) :: dt_tracers
real, allocatable, dimension(:,:,:) ::                                        &
     diff_m,               &   ! momentum diffusion coeff.
     diff_t,               &   ! temperature diffusion coeff.
     tdtlw,                &   ! place holder. appears in calling arguments of vert_turb_driver but not used unless do_edt=.true. -- pjp
     diss_heat,            &   ! heat dissipated by vertical diffusion
     non_diff_dt_ug,       &   ! zonal wind tendency except from vertical diffusion
     non_diff_dt_vg,       &   ! merid. wind tendency except from vertical diffusion
     non_diff_dt_tg,       &   ! temperature tendency except from vertical diffusion
     non_diff_dt_qg,       &   ! moisture tendency except from vertical diffusion
     conv_dt_tg,           &   ! temperature tendency from convection
     conv_dt_qg,           &   ! moisture tendency from convection
     cond_dt_tg,           &   ! temperature tendency from condensation
     cond_dt_qg                ! moisture tendency from condensation

logical, allocatable, dimension(:,:) ::                                       &
     avail,                &   ! generate surf. flux (all true)
     land,                 &   ! land points (all false)
     coldT,                &   ! should precipitation be snow at this point
     convect                   ! place holder. appears in calling arguments of vert_turb_driver but not used unless do_entrain=.true. -- pjp

real, allocatable, dimension(:,:) ::                                          &
     klzbs,                &   ! stored level of zero buoyancy values
     cape,                 &   ! convectively available potential energy
     cin,                  &   ! convective inhibition (this and the above are before the adjustment)
     invtau_q_relaxation,  &   ! temperature relaxation time scale
     invtau_t_relaxation,  &   ! humidity relaxation time scale
     rain,                 &   !
     snow

real, allocatable, dimension(:,:,:) :: &
     t_ref,          &   ! relaxation temperature for bettsmiller scheme
     q_ref               ! relaxation moisture for bettsmiller scheme

integer, allocatable, dimension(:,:) :: & convflag ! indicates which qe convection subroutines are used

integer ::           &
     id_diff_dt_ug,  &   ! zonal wind tendency from vertical diffusion
     id_diff_dt_vg,  &   ! merid. wind tendency from vertical diffusion
     id_diff_dt_tg,  &   ! temperature tendency from vertical diffusion
     id_diff_dt_qg,  &   ! moisture tendency from vertical diffusion
     id_conv_rain,   &   ! rain from convection
     id_cond_rain,   &   ! rain from condensation
     id_conv_dt_tg,  &   ! temperature tendency from convection
     id_conv_dt_qg,  &   ! temperature tendency from convection
     id_cond_dt_tg,  &   ! temperature tendency from convection
     id_cond_dt_qg       ! temperature tendency from convection

integer :: id_phalf, id_pfull, id_zhalf, id_zfull
integer :: id_drag, id_ustar, id_bstar
integer :: id_fluxt, id_ediff
integer :: id_bucket_depth         ! bucket depth vaiable for output added by LJJ
integer :: id_bucket_depth_lh, id_bucket_depth_cond, id_bucket_depth_conv
integer :: id_total_rain, id_rehum
integer :: id_total_cond, id_cld_evap, id_flux_vcld

real, allocatable, dimension(:,:,:) :: cld_dt_tg, cld_dt_qg !tendency from cloud evaporation
real, allocatable, dimension(:,:,:) :: rehum
real, allocatable, dimension(:,:) :: ediff
real, allocatable, dimension(:,:) :: bucket_depth, dt_bucket, &
     depth_change_lh,      &   ! tendency in bucket depth due to latent heat transfer     ! added by LJJ
     depth_change_cond,    &   ! tendency in bucket depth due to condensation rain        ! added by LJJ
     depth_change_conv         ! tendency in bucket depth due to convection rain ! added by LJJ
 

type(surf_diff_type) :: Tri_surf ! used by gcm_vert_diff
integer i,j, isc, iec, jsc, jec, num_levels

!=================================================================================================================================

contains

!#######################################################################

  subroutine atmosphere_init ( Time_init, Time, Time_step )

    type (time_type), intent(in) :: Time_step
    type (time_type), intent(in) :: Time_init
    type (time_type), intent(in) :: Time

    ! local:
    integer :: axes(5) !4)
    integer :: ss, ds
    integer :: ntiles=1
!    integer i,j, isc, iec, jsc, jec, num_levels
    integer :: unit, ierr, io
    real pp(2)


  !----- write version and namelist to log file -----

    call write_version_number ( version, tag )

  !---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

!FD, read namelist
unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=atmosphere_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'atmosphere_nml')
enddo
10 call close_file (unit)
if ( mpp_pe() == mpp_root_pe() )   write (stdlog(), nml=atmosphere_nml)


  !----- initialize FV dynamical core -----
!   cold_start = (.not.file_exist('INPUT/fv_rst.res.nc').and. .not.file_exist('INPUT/'//fms_tracers_file))
    cold_start = (.not.file_exist('INPUT/fv_core.res.nc'))

    allocate(Atm(ntiles))
    call fv_init(Atm(:),dt_atmos)  ! allocates Atm components

    Atm(1)%moist_phys = .false.

    ! Init model data
         call timing_on('fv_restart')
    call fv_restart(domain, Atm, dt_atmos, seconds, days, cold_start, grid_type)
         call timing_off('fv_restart')

#ifdef PERTURB_IC
    isc = Atm(1)%isc
    iec = Atm(1)%iec
    jsc = Atm(1)%jsc
    jec = Atm(1)%jec
#ifdef SW_DYNAMICS
! TEST CASE-7
!     if ( .not. cold_start ) then
         do j=jsc,jec
            do i=isc,iec
               pp(1) = Atm(1)%agrid(i,j,1)
               pp(2) = Atm(1)%agrid(i,j,2)
               Atm(1)%delp(i,j,1) = Atm(1)%delp(i,j,1) + 120.*grav*cos(pp(2)) *  &
               exp( -(3.*(pp(1)-pi))**2 ) * exp( -(15.*(pp(2)-pi/4.))**2 )
            enddo
         enddo
!     endif
#endif
#endif

! Need to implement time in fv_restart files
!   if ( .not. cold_start ) then 
! !
! ! Check consistency in Time
! !
!     fv_time = set_time (seconds, days)
!     call get_time (Time, ss,  ds)
!
!     if(seconds /= ss .or. days /= ds) then
!        unit = stdout()
!       write(unit,*) 'FMS:', ds, ss
!       write(unit,*) 'FV:', days, seconds
!       call error_mesg('FV_init:','Time inconsistent between fv_rst and INPUT/atmos_model.res',FATAL)
!     endif
!   else
      fv_time = time
!   endif



    call fv_diag_init(Atm, axes, Time, Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, p_ref)
    !call lin_cld_microphys_init(axes, Time)

!   if( nlev > 1 ) call hs_forcing_init ( axes, Time )

!-----------------------------------------------------------------------
!FD, ts_init
    isc = Atm(1)%isc
    iec = Atm(1)%iec
    jsc = Atm(1)%jsc
    jec = Atm(1)%jec
    num_levels = Atm(1)%npz

!FD, PBL
allocate(dt_ua(Atm(1)%isd:Atm(1)%ied, &
        Atm(1)%jsd:Atm(1)%jed, num_levels))
allocate(dt_va(Atm(1)%isd:Atm(1)%ied, &
        Atm(1)%jsd:Atm(1)%jed, num_levels))

allocate(tg_tmp      (isc:iec, jsc:jec, num_levels))
allocate(qg_tmp      (isc:iec, jsc:jec, num_levels))
allocate(dt_ug       (isc:iec, jsc:jec, num_levels))
allocate(dt_vg       (isc:iec, jsc:jec, num_levels))
allocate(dt_tg       (isc:iec, jsc:jec, num_levels))
allocate(dt_tracers  (isc:iec, jsc:jec, num_levels, Atm(1)%ncnst))
allocate(non_diff_dt_tg (isc:iec, jsc:jec, num_levels))
allocate(non_diff_dt_ug (isc:iec, jsc:jec, num_levels))

allocate(p_full      (isc:iec, jsc:jec, num_levels))
allocate(z_full      (isc:iec, jsc:jec, num_levels))
allocate(p_half      (isc:iec, jsc:jec, num_levels+1))
allocate(z_half      (isc:iec, jsc:jec, num_levels+1)); z_half = 0.
allocate(z_surf      (isc:iec, jsc:jec)); z_surf = 0.0
!allocate(t_surf      (isc:iec, jsc:jec))
allocate(q_surf      (isc:iec, jsc:jec)); q_surf = 0.0
allocate(u_surf      (isc:iec, jsc:jec)); u_surf = 0.0
allocate(v_surf      (isc:iec, jsc:jec)); v_surf = 0.0
allocate(rough_mom   (isc:iec, jsc:jec)); rough_mom = roughness_mom
allocate(rough_heat  (isc:iec, jsc:jec)); rough_heat = roughness_heat
allocate(rough_moist (isc:iec, jsc:jec)); rough_moist = roughness_moist
allocate(gust        (isc:iec, jsc:jec)); gust = 1.0
allocate(z_pbl       (isc:iec, jsc:jec))
allocate(flux_t      (isc:iec, jsc:jec))
allocate(flux_q      (isc:iec, jsc:jec))
allocate(flux_r      (isc:iec, jsc:jec))
allocate(flux_u      (isc:iec, jsc:jec))
allocate(flux_v      (isc:iec, jsc:jec))
allocate(drag_m      (isc:iec, jsc:jec))
allocate(drag_t      (isc:iec, jsc:jec))
allocate(drag_q      (isc:iec, jsc:jec))
allocate(w_atm       (isc:iec, jsc:jec))
allocate(ustar       (isc:iec, jsc:jec))
allocate(bstar       (isc:iec, jsc:jec))
allocate(qstar       (isc:iec, jsc:jec))
allocate(dhdt_surf   (isc:iec, jsc:jec))
allocate(dedt_surf   (isc:iec, jsc:jec))
allocate(dedq_surf   (isc:iec, jsc:jec))
allocate(drdt_surf   (isc:iec, jsc:jec))
allocate(dhdt_atm    (isc:iec, jsc:jec))
allocate(dedq_atm    (isc:iec, jsc:jec))
allocate(dtaudv_atm  (isc:iec, jsc:jec))
allocate(dtaudu_atm  (isc:iec, jsc:jec))
allocate(land        (isc:iec, jsc:jec)); land = .false.
allocate(avail       (isc:iec, jsc:jec)); avail = .true.
allocate(fracland    (isc:iec, jsc:jec)); fracland = 0.0
allocate(rough       (isc:iec, jsc:jec))
allocate(diff_t      (isc:iec, jsc:jec, num_levels))
allocate(diff_m      (isc:iec, jsc:jec, num_levels))
allocate(diss_heat   (isc:iec, jsc:jec, num_levels))

allocate(tdtlw       (isc:iec, jsc:jec, num_levels)); tdtlw = 0.0
allocate(convect     (isc:iec, jsc:jec)); convect = .false.
allocate(ediff       (isc:iec, jsc:jec)); ediff = 0.       

!moist
allocate(conv_dt_tg  (isc:iec, jsc:jec, num_levels))
allocate(conv_dt_qg  (isc:iec, jsc:jec, num_levels))
allocate(cond_dt_tg  (isc:iec, jsc:jec, num_levels))
allocate(cond_dt_qg  (isc:iec, jsc:jec, num_levels))
allocate(coldT        (isc:iec, jsc:jec)); coldT = .false.
allocate(klzbs        (isc:iec, jsc:jec))
allocate(cape         (isc:iec, jsc:jec))
allocate(cin          (isc:iec, jsc:jec))
allocate(invtau_q_relaxation  (isc:iec, jsc:jec))
allocate(invtau_t_relaxation  (isc:iec, jsc:jec))
allocate(rain         (isc:iec, jsc:jec)); rain = 0.0
allocate(snow         (isc:iec, jsc:jec)); snow = 0.0
allocate(convflag     (isc:iec, jsc:jec))
allocate(t_ref (isc:iec, jsc:jec, num_levels)); t_ref = 0.0
allocate(q_ref (isc:iec, jsc:jec, num_levels)); q_ref = 0.0
allocate(conv_cloud   (isc:iec, jsc:jec, num_levels)) !cloud mass fraction
allocate(rain_profile (isc:iec, jsc:jec, num_levels)) !rain mass fraction
allocate(cld_dt_tg    (isc:iec, jsc:jec, num_levels))
allocate(cld_dt_qg    (isc:iec, jsc:jec, num_levels))
allocate(rehum        (isc:iec, jsc:jec, num_levels))

!bucket model
allocate(dt_bucket        (isc:iec, jsc:jec)) 
allocate(bucket_depth     (isc:iec, jsc:jec))        ! added by LJJ
allocate(depth_change_lh  (isc:iec, jsc:jec)); depth_change_lh=0.0
allocate(depth_change_cond(isc:iec, jsc:jec)); depth_change_cond=0.0
allocate(depth_change_conv(isc:iec, jsc:jec)); depth_change_conv=0.0

!FD, ts_init
    if (cold_start) then
       do j=jsc,jec
          do i=isc,iec
             Atm(1)%ts(i,j) = Atm(1)%pt(i,j,Atm(1)%npz) + 1.
          enddo
       enddo
       Atm(1)%q(:,:,:,1) = 1e-6  !has to be zero for pure CO2 rad
       Atm(1)%q(:,:,:,4) = 1e-9  !used as the cloud mass fraction
!       if (bucket)  bucket_depth(:,:) = init_bucket_depth   
      
       Atm(1)%u_srf(:,:) = 0.
        if (bucket)  then
           bucket_depth(:,:) = init_bucket_depth
           bucket_depth(:,:) = Atm(1)%agrid(isc:iec,jsc:jec,2) *1.
!           bucket_depth(isc:iec,jsc:jec) = cos(Atm(1)%agrid(isc:iec,jsc:jec,2)) * &
!                cos(Atm(1)%agrid(isc:iec,jsc:jec,1) - 270. / 180.*pi)
!           !0.34~sin(20./180*pi)
!           !where (bucket_depth .gt. 0.94 .or. bucket_depth .lt. 0.0) 
!           where (bucket_depth .gt.     80./180.*pi) 
           where (bucket_depth .lt.    -80./180.*pi) 
             bucket_depth(:,:) = init_bucket_depth
           elsewhere
             bucket_depth(:,:) = 0.0
           endwhere
        endif

    else
!u_srf here is used for restart
!        Atm(1)%ts(:,:) = Atm(1)%u_srf(:,:)
!        if (read_bucket) call get_latlon_bucket(Atm)
!        if (bucket) bucket_depth(isc:iec, jsc:jec) = Atm(1)%v_srf(isc:iec,jsc:jec)

        !Atm(1)%q(:,:,:,1) = 4.5e-9
        if (bucket) then
           bucket_depth(isc:iec, jsc:jec) = Atm(1)%v_srf(isc:iec,jsc:jec) !q(isc:iec,jsc:jec,1,2) 
           dt_bucket   (isc:iec, jsc:jec) = Atm(1)%u_srf(isc:iec,jsc:jec) 
           Atm(1)%u_srf(isc:iec,jsc:jec)  = 0.

           if (bucket_iter_res) then
           bucket_i = g_sum(bucket_depth(isc:iec,jsc:jec), isc, iec, jsc, jec,Atm(1)%ng, area, 0)
           ! bucket time tendency within $days iteration, 
           !dt_bucket (365*100) /($days)
           bucket_depth(:,:) = bucket_depth(:,:) + dt_bucket(:,:) *(365./200) *bucket_scale
           where (bucket_depth <= 0.) bucket_depth = 0.
           where (bucket_depth >=40.) bucket_depth = 40.   
           !normalization
           bucket_f = g_sum(bucket_depth(isc:iec,jsc:jec), isc, iec, jsc, jec,Atm(1)%ng, area, 0)
           bucket_depth(:,:) = bucket_depth(:,:) / bucket_f * bucket_i
           endif
        endif

    endif

!LBL
do i=1,num_levels
!   p_half(:,:,i) = Atm(1)%pe(isc:iec, i, jsc:jec)
   p_full(:,:,i) = Atm(1)%delp(isc:iec,jsc:jec,i) / &
        (Atm(1)%peln(:,i+1,:) - Atm(1)%peln(:,i,:))
enddo

!Atm(1)%q(:,:,:,1) = 1e-6  !has to be zero for pure CO2 rad
!!for 1D comparison
!Atm(1)%ts(:,:) = 300.
!Atm(1)%pt(isc:iec, jsc:jec, :) = 300* &
!        (p_full(isc:iec, jsc:jec, :) / ps0)**kappa

call radiance_init(isc, iec, jsc, jec, &
                num_levels, axes, Time, &
                Atm(1)%agrid(isc:iec,jsc:jec,:), &
                Atm(1)%q(isc:iec,jsc:jec,:,1), &
                Atm(1)%pt(isc:iec,jsc:jec,:), & 
                Atm(1)%ps(isc:iec,jsc:jec), &
                p_full(:,:,:))
     !from 1 to npz

call mixed_layer_init(isc, iec, jsc, jec, num_levels, t_surf, axes(1:4), Time) 
call vert_diff_init (Tri_surf, iec-isc+1, jec-jsc+1, num_levels, .true., do_virtual)
call vert_turb_driver_init (Atm(1)%agrid(isc:iec,jsc:jec,1), Atm(1)%agrid(isc:iec,jsc:jec,2), & 
        iec-isc+1,jec-jsc+1, &
        num_levels,axes(1:4),Time, doing_edt, doing_entrain)
!isc:iec not boundary

id_phalf = register_diag_field(mod_name, 'ph',        &
     axes((/1,2,4/)), Time, 'half pressure','Pa')
id_pfull = register_diag_field(mod_name, 'pf',        &
     axes((/1,2,3/)), Time, 'full pressure','Pa')
id_zhalf = register_diag_field(mod_name, 'zhalf',        &
     axes((/1,2,4/)), Time, 'half height','m')
id_zfull = register_diag_field(mod_name, 'zfull',        &
     axes((/1,2,3/)), Time, 'full height','m')

id_fluxt = register_diag_field(mod_name, 'fluxt',        &
     axes((/1,2/)), Time, 'sensible heat','W/m**2')
id_ediff = register_diag_field(mod_name, 'ediff',        &
     axes((/1,2/)), Time, 'ediff','W/m**2')

id_drag = register_diag_field(mod_name, 'dragm',        &
     axes((/1,2/)), Time, 'dragm','W/m**2')
id_ustar = register_diag_field(mod_name, 'ustar',        &
     axes((/1,2/)), Time, 'ustar','W/m**2')
id_bstar = register_diag_field(mod_name, 'bstar',        &
     axes((/1,2/)), Time, 'bstar','W/m**2')

call lscale_cond_init()
id_cond_dt_qg = register_diag_field(mod_name, 'dt_qg_condensation',        &
     axes(1:3), Time, 'Moisture tendency from condensation','kg/kg/s')
id_cond_dt_tg = register_diag_field(mod_name, 'dt_tg_condensation',        &
     axes(1:3), Time, 'Temperature tendency from condensation','K/s')
id_cond_rain = register_diag_field(mod_name, 'condensation_rain',          &
     axes(1:2), Time, 'Rain from condensation','kg/m/m/s')

id_total_rain = register_diag_field(mod_name, 'total_rain',          &
     axes(1:2), Time, 'Rain rate','watt/m/m')
id_rehum      = register_diag_field(mod_name, 'rehum',          &
     axes(1:3), Time, 'Relative humidity','none')
id_total_cond = register_diag_field(mod_name, 'total_cond',          &
     axes(1:3), Time, 'Condensation rate','kg/kg/s')
id_cld_evap   = register_diag_field(mod_name, 'cld_evap',          &
     axes(1:3), Time, 'Evaporation rate of clouds','kg/kg/s')
id_flux_vcld  = register_diag_field(mod_name, 'flux_vcld',          &
     axes(1:3), Time, 'Meridional cloud mass flux','m/s kg/kg')

if (do_bettsmiller .or. do_manabe) then
   call qe_moist_convection_init()
   id_conv_dt_qg = register_diag_field(mod_name, 'dt_qg_convection',          &
        axes(1:3), Time, 'Moisture tendency from convection','kg/kg/s')
   id_conv_dt_tg = register_diag_field(mod_name, 'dt_tg_convection',          &
        axes(1:3), Time, 'Temperature tendency from convection','K/s')
   id_conv_rain = register_diag_field(mod_name, 'convection_rain',            &
        axes(1:2), Time, 'Rain from convection','kg/m/m/s')
endif

if(bucket) then
  id_bucket_depth = register_diag_field(mod_name, 'bucket_depth', &         ! added by LJJ
       axes(1:2), Time, 'Depth of surface reservoir', 'm')
  id_bucket_depth_conv = register_diag_field(mod_name, 'bucket_depth_conv', &         ! added by LJJ
        axes(1:2), Time, 'Tendency of bucket depth induced by Convection', 'm/s')
  id_bucket_depth_cond = register_diag_field(mod_name, 'bucket_depth_cond', &         ! added by LJJ
        axes(1:2), Time, 'Tendency of bucket depth induced by Condensation','m/s')
  id_bucket_depth_lh = register_diag_field(mod_name, 'bucket_depth_lh', &         ! added by LJJ
        axes(1:2), Time, 'Tendency of bucket depth induced by LH', 'm/s')

endif


  end subroutine atmosphere_init


!#######################################################################

  subroutine atmosphere (Time)
    type(time_type), intent(in) :: Time

    real:: zvir
    real:: time_total
    real:: tau_winds, tau_press, tau_temp

    integer i,j,k!, isc, iec, jsc, jec
    real :: hr_top

#ifdef NUDGE_IC
    tau_winds =  1. * 3600.
    tau_press = -1.
    tau_temp  = -1.
#else
    tau_winds = -1.
    tau_press = -1.
    tau_temp  = -1.
#endif

    fv_time = Time + Time_step_atmos
    call get_time (fv_time, seconds,  days)

    time_total = days*SECONDS_PER_DAY + seconds

    if ( tau_winds>0. .or. tau_press>0. .or. tau_temp>0. )     &
    call  fv_nudge(Atm(1)%npz, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, Atm(1)%ng, &
                   Atm(1)%u, Atm(1)%v, Atm(1)%delp, Atm(1)%pt, dt_atmos,    &
                   tau_winds, tau_press, tau_temp)

  !---- call fv dynamics -----
  !  if ( adiabatic .or. Atm(1)%do_Held_Suarez ) then
  !       zvir = 0.         ! no virtual effect
  !  else
         zvir = rvgas/rdgas - 1.
  !  endif

    call set_domain(Atm(1)%domain)  ! needed for diagnostic output done in fv_dynamics
    call timing_on('fv_dynamics')
    call fv_dynamics(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%ncnst-Atm(1)%pnats,         &
                     Atm(1)%ng, dt_atmos, Atm(1)%consv_te,                                  &
                     Atm(1)%fill, Atm(1)%reproduce_sum, kappa, cp_air, zvir,                &
                     Atm(1)%ks, Atm(1)%ncnst, Atm(1)%n_split, Atm(1)%q_split,               &
                     Atm(1)%u, Atm(1)%v, Atm(1)%um, Atm(1)%vm,                              &
                     Atm(1)%w, Atm(1)%delz, Atm(1)%hydrostatic,                             &
                     Atm(1)%pt, Atm(1)%delp, Atm(1)%q, Atm(1)%ps,                           &
                     Atm(1)%pe, Atm(1)%pk, Atm(1)%peln, Atm(1)%pkz, Atm(1)%phis,            &
                     Atm(1)%omga, Atm(1)%ua, Atm(1)%va, Atm(1)%uc, Atm(1)%vc,               &
                     Atm(1)%ak, Atm(1)%bk, Atm(1)%mfx, Atm(1)%mfy,                          &
                     Atm(1)%cx, Atm(1)%cy, Atm(1)%ze0, Atm(1)%hybrid_z, time_total)
    call timing_off('fv_dynamics')

!DF LBL===================================================================
isc = Atm(1)%isc
iec = Atm(1)%iec
jsc = Atm(1)%jsc
jec = Atm(1)%jec
num_levels = Atm(1)%npz
!p_atm = Atm(1)%delp(isc:iec,jsc:jec,num_levels) / &
!        (Atm(1)%peln(:,num_levels+1,:) - &
!        Atm(1)%peln(:,num_levels,:) )
!z_atm = (p_atm / Atm(1)%ps(isc:iec,jsc:jec)) * &
!        rdgas * Atm(1)%pt(isc:iec,jsc:jec,num_levels) /grav
do i=1,num_levels
   p_half(:,:,i) = Atm(1)%pe(isc:iec, i, jsc:jec)
   p_full(:,:,i) = Atm(1)%delp(isc:iec,jsc:jec,i) / &
        (Atm(1)%peln(:,i+1,:) - Atm(1)%peln(:,i,:))
enddo
p_half(:,:,num_levels+1) = Atm(1)%ps(isc:iec, jsc:jec)
z_half(:,:,1:num_levels) = Atm(1)%delp(isc:iec,jsc:jec,:) / &
        p_full(:,:,:) * rdgas /grav * (1+zvir*Atm(1)%q(isc:iec,jsc:jec,:,1)) *&
        Atm(1)%pt(isc:iec,jsc:jec,:)
z_full(:,:,:) = (p_half(:,:,2:num_levels+1) - p_full(:,:,:)) / &
        p_full(:,:,:) * rdgas /grav * (1+zvir*Atm(1)%q(isc:iec,jsc:jec,:,1)) *&
        Atm(1)%pt(isc:iec,jsc:jec,:)

do i=num_levels, 1, -1
   z_half(:,:,i) = z_half(:,:,i) + z_half(:,:,i+1)
   z_full(:,:,i) = z_full(:,:,i) + z_half(:,:,i+1)
enddo

!!DF DRY CONVEC=================================================================== 
!do i=isc, iec
!   do j=jsc, jec
!      call dry_convection(Atm(1)%pt(i,j,:), p_full(i,j,:), p_half(i,j,:), &
!      tg_tmp(i,j,:))
!   enddo
!enddo
!Atm(1)%pt(isc:iec, jsc:jec, :) = tg_tmp(isc:iec, jsc:jec, :)*1.


!if (mod(seconds, 14400) .eq. int(dt_atmos)) then !time_step_size
!       call radiance_update(Atm(1)%isc, Atm(1)%iec, &
!                Atm(1)%jsc, Atm(1)%jec, &
!                Atm(1)%ng,  Atm(1)%ncnst, Atm(1)%npz, &
!                fv_time, Atm(1)%agrid, &
!                Atm(1)%q, Atm(1)%pt, & 
!                Atm(1)%ts,Atm(1)%delp, Atm(1)%peln, &
!                non_diff_dt_tg       , Atm(1)%dTs_rad)
!                Atm(1)%dTrad, Atm(1)%dTs_rad)

!the 4th dimension of q is the cloud mass concentration
if (mod(seconds, deltat_rad) .eq. int(dt_atmos)) then !time_step_size
       call radiance_update(isc, iec, jsc, jec, &
                num_levels, fv_time, &
                Atm(1)%agrid(isc:iec,jsc:jec,:), &
                Atm(1)%q(isc:iec,jsc:jec,:,1), &
                Atm(1)%q(isc:iec,jsc:jec,:,4), &
                Atm(1)%pt(isc:iec,jsc:jec,:), & 
                Atm(1)%ts(isc:iec,jsc:jec), &
                bucket_depth(:,:), &
                p_half(:,:,:), p_full(:,:,:), &
                non_diff_dt_tg(:,:,:), &
                Atm(1)%dTs_rad(isc:iec,jsc:jec))

!!stablize the top layer, without breaking the energy conservation
!hr_top = globalsum(non_diff_dt_tg(isc:iec,jsc:jec,1), &
!        Atm(1)%npx, Atm(1)%npy, isc, iec, jsc, jec)
!non_diff_dt_tg(isc:iec, jsc:jec, 1) = hr_top
!if(gid==masterproc) write(6,*) 'TOA avergaed heating rate',hr_top 

endif

!update ta,ts
!    do i=Atm(1)%isc,Atm(1)%iec
!    do j=Atm(1)%jsc,Atm(1)%jec
!       Atm(1)%pt(i,j,:) = Atm(1)%pt(i,j,:) + &
!        Atm(1)%dTrad(i,j,:) * dt_atmos
!    enddo
!    enddo
!    Atm(1)%ts(:,:) = Atm(1)%ts(:,:) &
!        + Atm(1)%dTs_rad(:,:) * dt_atmos / 5.

!DF PBL===================================================================
dt_tg = non_diff_dt_tg *1. !radiative heating rate now
dt_ug = 0.; dt_vg = 0.; dt_tracers = 0.
dt_ua = 0.; dt_va = 0.
dt_bucket = 0.

call surface_flux(                                                          &
        Atm(1)%pt(isc:iec,jsc:jec,num_levels),                              &
       Atm(1)%q(isc:iec,jsc:jec,num_levels,1),                              &
        Atm(1)%ua(isc:iec,jsc:jec,num_levels),                              &
        Atm(1)%va(isc:iec,jsc:jec,num_levels),                              &
                       p_full(:,:,num_levels),                              &
           z_full(:,:,num_levels)-z_surf(:,:),                              &
                   Atm(1)%ps(isc:iec,jsc:jec),                              &
                   Atm(1)%ts(isc:iec,jsc:jec),                              &
                   Atm(1)%ts(isc:iec,jsc:jec),                              &
                                  q_surf(:,:),                              & ! is intent(inout)
                            bucket_depth(:,:),                              &
! added by LJJ
                         depth_change_lh(:,:),                              &
! added by LJJ
                       depth_change_conv(:,:),                              &
! added by LJJ
                       depth_change_cond(:,:),                              &
! added by LJJ
                                  u_surf(:,:),                              &
                                  v_surf(:,:),                              &
                               rough_mom(:,:),                              &
                              rough_heat(:,:),                              &
                             rough_moist(:,:),                              &
                               rough_mom(:,:),                              & ! using rough_mom in place of rough_scale -- pjp
                                    gust(:,:),                              &
                                  flux_t(:,:),                              & ! is intent(out)
                                  flux_q(:,:),                              & ! is intent(out)
                                  flux_r(:,:),                              & ! is intent(out)
                                  flux_u(:,:),                              & ! is intent(out)
                                  flux_v(:,:),                              & ! is intent(out)
                                  drag_m(:,:),                              & ! is intent(out)
                                  drag_t(:,:),                              & ! is intent(out)
                                  drag_q(:,:),                              & ! is intent(out)
                                   w_atm(:,:),                              & ! is intent(out)
                                   ustar(:,:),                              & ! is intent(out)
                                   bstar(:,:),                              & ! is intent(out)
                                   qstar(:,:),                              & ! is intent(out)
                               dhdt_surf(:,:),                              & ! is intent(out)
                               dedt_surf(:,:),                              & ! is intent(out)
                               dedq_surf(:,:),                              & ! is intent(out)
                               drdt_surf(:,:),                              & ! is intent(out)
                                dhdt_atm(:,:),                              & ! is intent(out)
                                dedq_atm(:,:),                              & ! is intent(out)
                              dtaudu_atm(:,:),                              & ! is intent(out)
                              dtaudv_atm(:,:),                              & ! is intent(out)
                                     dt_atmos,                              &
                                    land(:,:),                              &
                               .not.land(:,:),                              &
                                   avail(:,:)  )


  call vert_turb_driver(            1,                              1, &
                            Time,                 Time+Time_step_atmos, &
                              dt_atmos, tdtlw(:,:,:),    fracland(:,:), &
                 p_half(:,:,:),          p_full(:,:,:), &
                 z_half(:,:,:),          z_full(:,:,:), &
                            ustar(:,:),                     bstar(:,:), &
                            qstar(:,:),                     rough(:,:), &
          Atm(1)%agrid(isc:iec,jsc:jec,2),                convect(:,:), &
          Atm(1)%ua(isc:iec,jsc:jec,:),   Atm(1)%va(isc:iec,jsc:jec,:), &
          Atm(1)%pt(isc:iec,jsc:jec,:),                                 &
          Atm(1)%q(isc:iec,jsc:jec,:,1), Atm(1)%q(isc:iec,jsc:jec,:,:), &
          Atm(1)%ua(isc:iec,jsc:jec,:),   Atm(1)%va(isc:iec,jsc:jec,:), &
          Atm(1)%pt(isc:iec,jsc:jec,:),                                 &
          Atm(1)%q(isc:iec,jsc:jec,:,1), Atm(1)%q(isc:iec,jsc:jec,:,:), &
                          dt_ug(:,:,:),                   dt_vg(:,:,:), &
                          dt_tg(:,:,:),            dt_tracers(:,:,:,1), &
                   dt_tracers(:,:,:,:),                  diff_t(:,:,:), &
                         diff_m(:,:,:),                      gust(:,:), &
                            z_pbl(:,:) )
!
!! Don't zero these derivatives as the surface flux depends implicitly
!! on the lowest level values
!! However it should be noted that these derivatives do not take into
!! account the change in the Monin-Obukhov coefficients, and so are not
!! very accurate.
!
!!$   dtaudv_atm = 0.0
!!$   dhdt_atm   = 0.0
!!$   dedq_atm   = 0.0

!if (.false.) then

   call gcm_vert_diff_down (1, 1,  dt_atmos, &
             Atm(1)%ua(isc:iec,jsc:jec,:),   Atm(1)%va(isc:iec,jsc:jec,:), &
             Atm(1)%pt(isc:iec,jsc:jec,:),                                 &
             Atm(1)%q(isc:iec,jsc:jec,:,1), Atm(1)%q(isc:iec,jsc:jec,:,:), &
             diff_m(:,:,:),         diff_t(:,:,:),          p_half(:,:,:), &
             p_full(:,:,:),  z_full(:,:,:),     flux_u(:,:),  flux_v(:,:), &
                            dtaudu_atm(:,:),              dtaudv_atm(:,:), &
                            dt_ug(:,:,:),                    dt_vg(:,:,:), &
                            dt_tg(:,:,:),             dt_tracers(:,:,:,1), &
                            dt_tracers(:,:,:,:),         diss_heat(:,:,:), &
                            Tri_surf)

!
! update surface temperature
!
   call mixed_layer(                                                       &
                              Time,                                        &
                           Atm(1)%ts(:,:),                                 & ! t_surf is intent(inout)
                              flux_t(:,:),                                 &
                              flux_q(:,:),                                 &
                            fracland(:,:),                                 &
                                 dt_atmos,                                 &
                      Atm(1)%dTs_rad(:,:),                                 &
                            fracland(:,:),                                 & !dTs_rad = swd + lwn
                            Tri_surf,                                      & ! Tri_surf is intent(inout)
                           dhdt_surf(:,:),                                 &
                           dedt_surf(:,:),                                 &
                           dedq_surf(:,:),                                 &
                           drdt_surf(:,:),                                 &
                            dhdt_atm(:,:),                                 &
                            dedq_atm(:,:))

   call gcm_vert_diff_up (1, 1, dt_atmos, Tri_surf, dt_tg(:,:,:), dt_tracers(:,:,:,1), dt_tracers(:,:,:,:))

!endif

!DF Check Energy===================================================================
!non_diff_dt_ug = ((dt_tg - non_diff_dt_tg) * cp_air + &
!        Atm(1)%ua(isc:iec,jsc:jec,:) * dt_ug + &
!        Atm(1)%va(isc:iec,jsc:jec,:) * dt_vg) *&
!        Atm(1)%delp(isc:iec,jsc:jec,:) /grav
!ediff = 0.
!do i=1,num_levels
!   ediff = ediff + non_diff_dt_ug(:,:,i)
!enddo
!DF Betts-Miller convection ===================================================
if (do_bettsmiller) then
   rain = 0.0; snow = 0.0
   call qe_moist_convection (dt_atmos,    Atm(1)%pt(isc:iec,jsc:jec,:),      &
        Atm(1)%q(isc:iec,jsc:jec,:,1),                     p_full     ,      &
                          p_half     ,                           coldT,      &
                                 rain,                            snow,      &
                           conv_dt_tg,                      conv_dt_qg,      &
                                q_ref,                        convflag,      &
                                klzbs,                            cape,      &
                                  cin,             invtau_q_relaxation,      &
                  invtau_t_relaxation,                           t_ref)
   tg_tmp = conv_dt_tg + Atm(1)%pt(isc:iec,jsc:jec,:) 
   qg_tmp = conv_dt_qg + Atm(1)%q(isc:iec,jsc:jec,:,1) 
!  note the delta's are returned rather than the time derivatives
   conv_dt_tg = conv_dt_tg/dt_atmos
   conv_dt_qg = conv_dt_qg/dt_atmos
   rain       = rain/dt_atmos
   dt_tg      = dt_tg + conv_dt_tg
   dt_tracers(:,:,:,1) = dt_tracers(:,:,:,1) + conv_dt_qg
   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain > 0) used = send_data(id_conv_rain, rain, Time)
else if (do_manabe) then
   call moist_convection(Atm(1)%pt(isc:iec,jsc:jec,:), Atm(1)%q(isc:iec,jsc:jec,:,1), &
         p_full, p_half, &
         rain, t_ref, q_ref, conv_cloud) 
   tg_tmp     = t_ref *1.
   qg_tmp     = q_ref *1.
   conv_dt_tg = (t_ref - Atm(1)%pt(isc:iec,jsc:jec,:))/dt_atmos
   conv_dt_qg = (q_ref - Atm(1)%q(isc:iec,jsc:jec,:,1))/dt_atmos
   depth_change_conv = rain/dens_h2o
   rain       = rain / dt_atmos
   dt_tg      = dt_tg + conv_dt_tg
   dt_tracers(:,:,:,1) = dt_tracers(:,:,:,1) + conv_dt_qg
!   dt_tracers(:,:,:,4) = dt_tracers(:,:,:,4) - conv_dt_qg !conv_dt_qg <0
   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain > 0) used = send_data(id_conv_rain, rain*hlv, Time)
else 
   tg_tmp = Atm(1)%pt(isc:iec,jsc:jec,:) *1.
   qg_tmp = Atm(1)%q(isc:iec,jsc:jec,:,1) *1.
endif

!DF Lagescale condensation, t, q profiles after convection ===================================================
rain = 0.0
call lscale_cond (         tg_tmp,                   qg_tmp,               &
                    p_full(:,:,:),                   p_half(:,:,:),        &
                            coldT,                            rain,        &
                             snow,                      cond_dt_tg,        &
                       cond_dt_qg )                                
tg_tmp = cond_dt_tg + tg_tmp
qg_tmp = cond_dt_qg + qg_tmp
cond_dt_tg = cond_dt_tg/dt_atmos
cond_dt_qg = cond_dt_qg/dt_atmos
depth_change_cond = rain/dens_h2o  
rain       = rain/dt_atmos      
dt_tg      = dt_tg + cond_dt_tg
dt_tracers(:,:,:,1) = dt_tracers(:,:,:,1) + cond_dt_qg
!dt_tracers(:,:,:,4) = dt_tracers(:,:,:,4) - cond_dt_qg !cond_dt_qg <0
if(id_cond_dt_qg > 0) used = send_data(id_cond_dt_qg, cond_dt_qg, Time)
if(id_cond_dt_tg > 0) used = send_data(id_cond_dt_tg, cond_dt_tg, Time)
if(id_cond_rain  > 0) used = send_data(id_cond_rain, rain*hlv, Time)

!DF cloud and total_rain update=================================================================
!Atm(1)%q(isc:iec,jsc:jec,:,4) = Atm(1)%q(isc:iec,jsc:jec,:,4) + &
!        conv_cloud -cond_dt_qg !condensation rate -> cloud
rain = 0. 

rain_profile = Atm(1)%q(isc:iec,jsc:jec,:,4) / cloud_time !(3600)
Atm(1)%q(isc:iec,jsc:jec,:,4) = Atm(1)%q(isc:iec,jsc:jec,:,4)  &
        + conv_cloud -cond_dt_qg*dt_atmos & !condensation rate -> cloud
        - rain_profile *dt_atmos !autoconversion to precip
rain = sum(rain_profile(:,:,:)/grav * &
        (p_half(:,:,2:num_levels+1) - p_half(:,:,1:num_levels) ), 3)

!cloud evaporation, note that cld_dt is not tendency, but delta now
cld_dt_tg=0.; cld_dt_qg = 0.
call calculate_rehum(tg_tmp, qg_tmp, p_full, rehum)
call cloud_evap(dt_atmos, Atm(1)%q(isc:iec,jsc:jec,:,4), tg_tmp, qg_tmp, p_full, &
        rehum, cld_dt_tg, cld_dt_qg)
Atm(1)%q(isc:iec,jsc:jec,:,4) = Atm(1)%q(isc:iec,jsc:jec,:,4) - cld_dt_qg
dt_tg      = dt_tg + cld_dt_tg/dt_atmos
dt_tracers(:,:,:,1) = dt_tracers(:,:,:,1) + cld_dt_qg/dt_atmos

where (Atm(1)%q(isc:iec,jsc:jec,:,4) < 0.) Atm(1)%q(isc:iec,jsc:jec,:,4)=1e-20
if (any(Atm(1)%pt(isc:iec,jsc:jec,:) < 0.)) call error_mesg('cloud_evap:','neg T detected',FATAL)

!do i=isc,iec
!   do j=jsc,jec
!      do k=1,num_levels
!         if (Atm(1)%q(i,j,k,4) > 1e-6) then
!            rain(i,j) = rain(i,j) + (Atm(1)%q(i,j,k,4) - 1e-6)*&
!                (p_half(i,j,k+1)-p_half(i,j,k))/grav /dt_atmos
!            Atm(1)%q(i,j,k,4) = 1e-6
!         endif
!      enddo
!   enddo
!enddo
if(id_total_rain  >0) used = send_data(id_total_rain, rain*hlv, Time)
if(id_rehum       >0) used = send_data(id_rehum, rehum, Time)
if(id_total_cond  >0) used = send_data(id_total_cond, &
        0.-conv_dt_qg-cond_dt_qg, Time)
if(id_cld_evap    >0) used = send_data(id_cld_evap, cld_dt_qg/dt_atmos, Time)
if(id_flux_vcld   >0) used = send_data(id_flux_vcld, &
        Atm(1)%va(isc:iec,jsc:jec,:) *Atm(1)%q(isc:iec,jsc:jec,:,4), Time)
!DF UPDATE===================================================================
Atm(1)%ua(isc:iec,jsc:jec,:) = Atm(1)%ua(isc:iec,jsc:jec,:) + &
        dt_ug(:,:,:) * dt_atmos
Atm(1)%va(isc:iec,jsc:jec,:) = Atm(1)%va(isc:iec,jsc:jec,:) + &
        dt_vg(:,:,:) * dt_atmos
Atm(1)%pt(isc:iec,jsc:jec,:) = Atm(1)%pt(isc:iec,jsc:jec,:) + &
        dt_tg(:,:,:) * dt_atmos
Atm(1)%q(isc:iec,jsc:jec,:,1) = Atm(1)%q(isc:iec,jsc:jec,:,1) + &
        dt_tracers(:,:,:,1) * dt_atmos

!Atm(1)%q(isc:iec,jsc:jec,:,1) = 1e-6
!not now

                                                    call timing_on(' Update_dwinds')
    dt_ua(isc:iec,jsc:jec,:) = dt_ug
    dt_va(isc:iec,jsc:jec,:) = dt_vg
    call update_dwinds_phys(isc, iec, jsc, jec, &
        Atm(1)%isd, Atm(1)%ied, Atm(1)%jsd, Atm(1)%jed, &
        dt_atmos, dt_ua, dt_va, Atm(1)%u, Atm(1)%v)

                                                    call timing_off(' Update_dwinds')

!update the bucket depth
! added by LJJ
if(bucket) then
  
!   bucket_i = g_sum(bucket_depth(isc:iec,jsc:jec), isc, iec, jsc, jec, Atm(1)%ng, area, 0)

   ! bucket time tendency
   dt_bucket = depth_change_cond + depth_change_conv - depth_change_lh
   if (bucket_iter) bucket_depth(:,:) = bucket_depth(:,:) + dt_bucket 
   where (bucket_depth <= 0.) bucket_depth = 0.

!normalization
!   if (bucket_scale > 1) then
!      where (dt_bucket<0) bucket_depth = 0.
!      bucket_f = g_sum(bucket_depth(isc:iec,jsc:jec), isc, iec, jsc, jec, Atm(1)%ng, area, 0)
!      bucket_depth(:,:) = bucket_depth(:,:) / bucket_f * bucket_i
!      where (dt_bucket<0) bucket_depth = depth_change_cond + depth_change_conv 
!   endif

   if(id_bucket_depth > 0) used = send_data(id_bucket_depth,bucket_depth(:,:), Time)
   if(id_bucket_depth_conv > 0) used = send_data(id_bucket_depth_conv,depth_change_conv(:,:), Time)
   if(id_bucket_depth_cond > 0) used = send_data(id_bucket_depth_cond,depth_change_cond(:,:), Time)
   if(id_bucket_depth_lh > 0) used = send_data(id_bucket_depth_lh, dt_bucket, Time) !depth_change_lh(:,:), Time)
endif


if (id_fluxt > 0) used = send_data(id_fluxt, flux_t, Time)
if (id_ediff > 0) used = send_data(id_ediff, ediff, Time)

if (id_phalf > 0) used = send_data(id_phalf, p_half, Time)
if (id_pfull > 0) used = send_data(id_pfull, p_full, Time)
if (id_zhalf > 0) used = send_data(id_zhalf, z_half, Time)
if (id_zfull > 0) used = send_data(id_zfull, z_full, Time)

if (id_drag > 0) used = send_data(id_drag, drag_m, Time)
if (id_ustar > 0) used = send_data(id_ustar, ustar, Time)
if (id_bstar > 0) used = send_data(id_bstar, bstar, Time)

!DF UPDATE===================================================================
!    if(Atm(1)%npz /=1 .and. .not. adiabatic)then

!                                                         call timing_on('FV_PHYS')
!       call fv_phys(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%isc, Atm(1)%iec,          &
!                    Atm(1)%jsc, Atm(1)%jec, Atm(1)%ng, Atm(1)%ncnst,                     &
!                    Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%pt, Atm(1)%q, Atm(1)%pe,        &
!                    Atm(1)%delp, Atm(1)%peln, Atm(1)%pkz, dt_atmos,                      &
!                    Atm(1)%ua, Atm(1)%va, Atm(1)%phis, Atm(1)%agrid,                     &
!                    Atm(1)%ak, Atm(1)%bk, Atm(1)%ks, Atm(1)%ps, Atm(1)%pk,               &
!                    Atm(1)%u_srf, Atm(1)%v_srf, Atm(1)%ts, Atm(1)%delz,                  &
!                    Atm(1)%hydrostatic, Atm(1)%phys_hydrostatic, Atm(1)%oro, .true.,     &
!                    .false., p_ref, Atm(1)%fv_sg_adj, (mpp_pe()==mpp_root_pe()),         &
!                    Atm(1)%do_Held_Suarez, fv_time, time_total)
!                                                        call timing_off('FV_PHYS')
!    endif

    call nullify_domain()

  !---- diagnostics for FV dynamics -----

    call timing_on('FV_DIAG')

    call fv_diag(Atm, zvir, fv_time, Atm(1)%print_freq)

    call timing_off('FV_DIAG')


! used to restart Ts
!    Atm(1)%u_srf(:,:) = Atm(1)%ts(:,:)
!    if (bucket) Atm(1)%v_srf(isc:iec, jsc:jec) = bucket_depth(isc:iec, jsc:jec)
    if (bucket) then !v_srf: bucket_depth; u_srf: dt_bucket
!       Atm(1)%q(isc:iec,jsc:jec,1,2) = bucket_depth(isc:iec, jsc:jec)
       Atm(1)%v_srf(isc:iec,jsc:jec) = bucket_depth(isc:iec, jsc:jec)
       Atm(1)%u_srf(isc:iec,jsc:jec) = Atm(1)%u_srf(isc:iec,jsc:jec) + dt_bucket   (isc:iec, jsc:jec)
    endif

 end subroutine atmosphere


 subroutine atmosphere_end

    call get_time (fv_time, seconds,  days)

    !if ( Atm(1)%nwat==6 )    & 
    !call lin_cld_microphys_end
    
    if (do_bettsmiller) call qe_moist_convection_end

    call fv_end(Atm)
    deallocate(Atm)

  end subroutine atmosphere_end

 subroutine atmosphere_domain ( fv_domain )
 type(domain2d), intent(out) :: fv_domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos
        
   fv_domain = domain
        
 end subroutine atmosphere_domain

subroutine get_latlon_bucket(Atm)
type(fv_atmos_type), intent(inout):: Atm(:)
character(len=128) :: fname
real, allocatable:: lat(:), lon(:)
real, allocatable:: rdlat(:), rdlon(:)
real, allocatable:: bp0(:,:)
integer :: i, j, im, jm
integer tsize(4)
logical found
real:: a1, b1, c1, c2, c3, c4
real:: bpc
integer i1, i2, jc, i0, j0, iq
integer is, ie, js, je

! Read in lat-lon file
fname = Atm(1)%res_latlon_dynamics
call field_size(fname, 'temp', tsize, field_found=found)
if(gid==masterproc) write(6,*) 'Reconstruct cubed-sphere restart bucket depth', fname
!if(gid==masterproc) write(6,*) 'print_freq', Atm(1)%print_freq
          if ( found ) then
               im = tsize(1); jm = tsize(2)
               if(gid==masterproc)  write(6,*) 'External IC dimensions:', tsize
          endif
! Define the lat-lon coordinate:
          allocate (  lon(im) )
          allocate (  lat(jm) )
          allocate ( bp0(im,jm))
          allocate (rdlon(im) )
          allocate (rdlat(jm) )
          call read_data (fname, 'grid_yt', lat)
          call read_data (fname, 'grid_xt', lon)
          call read_data (fname, 'bucket_depth', bp0)
          do i=1,im
             lon(i) = lon(i) * (pi/180.)  ! lon(1) = 0.
          enddo
          do j=1,jm
             lat(j) = lat(j) * (pi/180.)
          enddo
  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))
  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo
! * Interpolate to cubed sphere cell center
is = Atm(1)%isc
ie = Atm(1)%iec
js = Atm(1)%jsc
je = Atm(1)%jec
do j=js,je
   do i=is,ie
       if ( Atm(1)%agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (Atm(1)%agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( Atm(1)%agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (Atm(1)%agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( Atm(1)%agrid(i,j,1)>=lon(i0) .and. Atm(1)%agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (Atm(1)%agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue
       if ( Atm(1)%agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( Atm(1)%agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( Atm(1)%agrid(i,j,2)>=lat(j0) .and. Atm(1)%agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (Atm(1)%agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue
       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) i,j,a1, b1
       endif
       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1
! Interpolated surface bucket depth
       bpc = c1*bp0(i1,jc  ) + c2*bp0(i2,jc  ) +    &
             c3*bp0(i2,jc+1) + c4*bp0(i1,jc+1)
       Atm(1)%v_srf(i,j) = bpc

   enddo
enddo
          deallocate (  lon )
          deallocate (  lat )
          deallocate ( bp0)
          deallocate (rdlon )
          deallocate (rdlat )
end subroutine get_latlon_bucket


end module atmosphere_mod



module mult_scat_core

  !--------------------------------------------------------
  !  
  ! Set up core multiple scattering routines.
  ! Code is two-stream only for now.
  !
  ! Robin Wordsworth (2018-2019)
  !
  !--------------------------------------------------------

  use constants_mod, only : pi

  implicit none
  private

  public :: setup_mult_scat_core, end_mult_scat_core
  public :: backscatter, angular_backscatter, backscatter_table
  public :: get_Iup_Idn


  !------- program options --------------  
  logical, parameter :: verbose = .false.

  save

contains

  ! ==================================================================================
  subroutine backscatter(g, b)
  implicit none
  real(8),intent(in) :: g
  real(8),intent(out):: b
  integer, parameter :: n=16
  integer m,i
  real(8) f, m0

  f  = g**(2*n)
  b  = 0.5 - 0.5*3*(g-f)/(1-f) /4.
  m0 = 1.
  do m=3,2*n-1, 2
     m0=m0*m/(m-1)
     b = b-0.5*(2*m+1) * (g**m-f) /(1-f)  &
        * ( (-1)**((m-1)/2.)/m/(m+1)*m0 )**2
  enddo
  return

  end subroutine backscatter

  ! ==================================================================================
  subroutine angular_backscatter(mu, g, b0)
  implicit none
  real(8),intent(in) :: mu, g
  real(8),intent(out):: b0
  integer, parameter :: n=16
  integer m,i
  real(8) b0_temp(0:2*n-1), f, m0

  call p_polynomial_value(1,2*n-1,mu,b0_temp)

  f  = g**(2*n)
  b0 = 0.5 - 0.5*3*(g-f)/(1-f) *b0_temp(1) /2.
  m0 = 1.
  do m=3,2*n-1, 2
     m0=m0*m/(m-1)
     b0 = b0-0.5*(2*m+1) * (g**m-f) /(1-f) *b0_temp(m) &
        * (-1)**((m-1)/2.)/m/(m+1)*m0
  enddo
  return

  end subroutine angular_backscatter

  ! ==================================================================================
  subroutine backscatter_table(g, mu0, b0)
  implicit none
  real(8),intent(in) :: g
  real(8),dimension(20), intent(out):: mu0, b0
  integer, parameter :: np=20
  integer, parameter :: n=16
  integer m,i
  real(8) b0_temp(np, 0:2*n-1), f, m0

  do i=1,np
     mu0(i)=-1+(i-1)*2./(np-1)
  enddo

  call p_polynomial_value(np,2*n-1,mu0,b0_temp)

  f  = g**(2*n)
  b0 = 0.5 - 0.5*3*(g-f)/(1-f) *b0_temp(:,1) /2.
  m0 = 1.
  do m=3,2*n-1, 2
     m0=m0*m/(m-1)
     b0 = b0-0.5*(2*m+1) * (g**m-f) /(1-f) *b0_temp(:,m) &
        * (-1)**((m-1)/2.)/m/(m+1)*m0
  enddo
  return

  end subroutine backscatter_table
  ! ==================================================================================

  subroutine setup_mult_scat_core()

    implicit none

    ! allocate all variables that require it here

    return

  end subroutine setup_mult_scat_core

  ! ==================================================================================

  subroutine end_mult_scat_core()

    implicit none

    return

  end subroutine end_mult_scat_core

  ! ==================================================================================

  subroutine Gam_fn(a,b,muav,Gam)

    ! See Thomas & Stamnes, S. 7.5.3, p. 258

    implicit none

    real(8), intent(in)  :: a    ! single scattering albedo []
    real(8), intent(in)  :: b    ! backscatter parameter []
    real(8), intent(in)  :: muav ! mean propagation angle cosine []
    real(8), intent(out) :: Gam  ! Gamma parameter []

    Gam = sqrt((1.0d0 - a)*(1.0d0 - a + 2.0d0*a*b))/muav

    return
  end subroutine Gam_fn

  ! ==================================================================================

  subroutine r_fn(a,b,r)

    ! See Thomas & Stamnes, S. 7.5.3, p. 258

    implicit none

    real(8), intent(in)  :: a    ! single scattering albedo []
    real(8), intent(in)  :: b    ! backscatter parameter []
    real(8), intent(out) :: r    ! albedo of semi-infinite homogenous slab []

    real(8) Num, Denom

    Num   = sqrt(1.0d0 - a + 2.0d0*a*b) - sqrt(1.0d0 - a)
    Denom = sqrt(1.0d0 - a + 2.0d0*a*b) + sqrt(1.0d0 - a)
    r     = Num/Denom

    return
  end subroutine r_fn

  ! ==================================================================================

  subroutine X_Y_fn(a,b,b0,muav,mu0,F0,X,Y)

    ! calculate the source terms from the direct beam to the multiple scattering part

    implicit none

    real(8), intent(in)  :: a    ! single scattering albedo []
    real(8), intent(in)  :: b    ! backscatter parameter []
    real(8), intent(in)  :: b0   ! backscatter parameter []
    real(8), intent(in)  :: muav ! mean propagation angle cosine []
    real(8), intent(in)  :: mu0  ! stellar zenith angle cosine []
    real(8), intent(in)  :: F0   ! incoming stellar flux [W/m2/cm^-1]
    real(8), intent(out) :: X    ! X parameter (source term for upward beam) [W/m2/cm^-1/sr]
    real(8), intent(out) :: Y    ! Y parameter (source term for downward beam) [W/m2/cm^-1/sr]

    real(8) Zc, Xplus, Xminus

!    Zc = a*F0*mu0 / ( 4.0d0*pi*muav**2.0d0*(1.0d0 - Gam**2.0d0*mu0**2.0d0) )
!    X  = +Zc*(muav - mu0)
!    Y  = -Zc*(muav + mu0)

        !b  = Gam ! backscatter parameter []
        !call angular_backscatter(mu0,0.9d0, b0)
        !print *, b0, b
        !b0 = 0.10
        !b0 = 0.5     !0.5d0*(1-0.9/1.9 *3*muav*mu0) ! backscatter parameter []
        Xplus  = a*F0/2/pi*b0
        Xminus = a*F0/2/pi*(1-b0)
        Zc = (1-a)*(1-a+2*a*b)-(muav/mu0)**2
        X  = (a*b*Xminus + (1-a+a*b-muav/mu0)*Xplus) /Zc
        Y  = (a*b*Xplus  + (1-a+a*b+muav/mu0)*Xminus)/Zc

    return

  end subroutine X_Y_fn

  ! ==================================================================================

  subroutine get_Iup_Idn(aL,bL,b0L,tauL,F0,As,muav,mu0,Iup,Idn,debug,nLay,nLev)

    implicit none

    real(8), intent(in)  :: aL(:) !(nLay)   ! single scattering albedo []
    real(8), intent(in)  :: bL(:) !(nLay)   ! backscatter parameter []
    real(8), intent(in)  :: b0L(:) !(nLay)   ! angular backscatter parameter []
    real(8), intent(in)  :: tauL(:) !(nLev) ! optical depth []
    real(8), intent(in)  :: F0         ! stellar flux [W/m2/cm^-1]
    real(8), intent(in)  :: As         ! surface albedo / Lambertian BDRF []
    real(8), intent(in)  :: mu0        ! stellar zenith angle cosine []
    real(8), intent(in)  :: muav       ! average propagation angle cosine []
    real(8), intent(out) :: Iup(:) !(nLev)  ! two-stream upwards irradiance [W/m2/cm^-1/sr]
    real(8), intent(out) :: Idn(:) !(nLev)  ! two-stream downwards irradiance [W/m2/cm^-1/sr]
    logical, intent(in)  :: debug
    integer, intent(in)  :: nLay, nLev

    integer iLayA, iLayB, j, k

  !------- general variables ----------
  integer :: N      ! total number of 2-stream elements to be solved []

  integer iLay, iLev                    ! do loop variables

  real(8) rL(nLay)                      ! semi-infinite slab albedo []
  real(8) sL(nLay)                      ! a function of rL []
  real(8) GamL(nLay)                    ! Gamma parameter []
  real(8) betaL(nLev)                   ! beta parameter  = Gamma*tau []

  real(8) X(nLay), Y(nLay)              ! up/down source coefficients [W/m2/sr CHECK]
  real(8) ExMu0(nLev)                   ! decay term for direct beam []
  real(8) Rdir                          ! reflected direct beam fraction []

    real(8) cL(nLay) ! 
    real(8) dL(nLay) ! 
    real(8) gL(nLay) ! 
    real(8) zL(nLay) ! 
    real(8) wL(nLay) ! 

    real(8),allocatable,dimension(:) :: aT(:) ! tri-di matrix 1st coefficient vector []
    real(8),allocatable,dimension(:) :: bT(:) ! tri-di matrix 2nd coefficient vector []
    real(8),allocatable,dimension(:) :: cT(:) ! tri-di matrix 3rd coefficient vector []
    real(8),allocatable,dimension(:) :: dT(:) ! tri-di matrix external vector []
    real(8),allocatable,dimension(:) :: xT(:) ! tri-di matrix solution vector []

    N = 2*nLay 
    allocate(aT (N))
    allocate(bT (N))
    allocate(cT (N))
    allocate(dT (N))
    allocate(xT (N))

    cL(:) = 0.0d0; dL(:) = 0.0d0; gL(:) = 0.0d0; zL(:) = 0.0d0; wL(:) = 0.0d0
    aT(:) = 0.0d0; bT(:) = 0.0d0; cT(:) = 0.0d0; dT(:) = 0.0d0; xT(:) = 0.0d0

    ! get Gamma, r and s 
    do iLay = 1,nLay
       call r_fn(aL(iLay),bL(iLay),rL(iLay))
       call Gam_fn(aL(iLay),bL(iLay),muav,GamL(iLay))
    end do
    sL(:) = rL(:) - 1.0d0/rL(:)

    ! get reflected direct beam fraction Rdir
    ExMu0(:) = exp(-tauL/mu0)
    Rdir     = As*mu0*F0*ExMu0(nLev)/(2*pi*muav) 
    !Rdir     = As*mu0*F0*ExMu0(nLev)/(pi*muav) 
    ! to be confirmed analytically that it's pi not 2*pi here (it works)...

    ! get the source functions
    do iLay = 1,nLay
       call X_Y_fn(aL(iLay),bL(iLay),b0L(iLay),muav,mu0,F0,X(iLay),Y(iLay))
       gL(iLay) = exp(-GamL(ilay)*(tauL(iLay+1)-tauL(iLay)) )
    end do

    ! get beta and its +ve/-ve exponents
    !betaL(1)      = GamL(1)*tauL(1)
    !betaL(2:nLev) = GamL(:)*tauL(2:nLev)
    ! Note: we are cutting a corner here, potentially, because GamL is defined 
    ! on the layers whereas tau is defined on the levels. Another choice would
    ! be to make GamL in betaL take the same layer as rL in each row of AMat. I
    ! will return to this issue later.
    !gL(:) = exp(-betaL(2:nLev) + betaL(1:nLev-1))

    ! calculate boundary condition / source term vector
    zL(1) = -Y(1)
    do j=2,(N - 1)
       iLev      = int(floor(real(j)/2)+1)
       iLayA     = iLev
       iLayB     = iLev-1
       zL(iLayA) = +(Y(iLayA) - Y(iLayB))*ExMu0(iLev) ! idn
       wL(iLayB) = +(X(iLayA) - X(iLayB))*ExMu0(iLev) ! iup
    end do
    wL(nLay) = (X(nLay) - As*Y(nLay))*ExMu0(nLev) - Rdir

    ! now solve matrix equation the tri-di way

    ! calculate matrix diagonals and the vector
    aT(1) = 0.0d0 
    bT(1) = 1.0d0
    cT(1) = rL(1) *gL(1)
    dT(1) = zL(1)
    do j=2,(N - 1)
       iLev = int(floor(real(j)/2.0d0)+1.0d0)
       iLay = iLev
       if(mod(j,2).eq.0)then ! is even
          !aT(j) = gL(iLev-1)
          !bT(j) = 0.0d0
          !cT(j) = -1.0d0
          !dT(j) = (wL(iLay-1) - zL(iLay)/rL(iLay))/sL(iLay)
          aT(j) = gL(iLev-1) *(rL(iLay-1) - 1./rL(iLay))
          bT(j) = (1.-rL(iLay-1)/rL(iLay)) 
          cT(j) = -sL(iLay)
          dT(j) = (wL(iLay-1) - zL(iLay)/rL(iLay)) !/sL(iLay)
       else ! is odd
          !aT(j) = 1.0d0 
          !bT(j) = 0.0d0
          !cT(j) = -gL(iLev)
          !dT(j) = (zL(iLay) - wL(iLay-1)/rL(iLay))/sL(iLay)
          aT(j) = sL(iLay-1) 
          bT(j) = -(1.-rL(iLay)/rL(iLay-1)) 
          cT(j) = -gL(iLev) *(rL(iLay) - 1./rL(iLay-1))
          dT(j) = (zL(iLay) - wL(iLay-1)/rL(iLay-1)) !/sL(iLay)
       end if
    end do
    aT(N) = (As - rL(nLay))*gL(nLay)
    bT(N) = As*rL(nLay) - 1.0d0
    cT(N) = 0.0d0
    dT(N) = wL(nLay)

    ! apply some voodoo numerics
    bT    = bT + 1.0d-16

    ! solve the problem using the Thomas algorithm
    call thomas_solve(N,aT,bT,cT,dT,xT)

    ! derive cL and dL from the big vector xT
    k = 1
    do j=1,N
       if(mod(j,2).eq.0)then ! is even
          dL(k) = xT(j)
          k = k + 1
       else ! is odd
          cL(k) = xT(j)
       end if
    end do

    ! get half-range upward irradiance [W/m2/sr]
    Iup(1:nLev-1) = cL(:)*rL(:)       + dL(:)*gL(:)              + X(:)*ExMu0(1:nLev-1)
    Iup(nLev)     = cL(nLay)*rL(nLay)*gL(nLay) +dL(nLay)         + X(nLay)*ExMu0(nLev)

    ! get half-range downward irradiance [W/m2/sr]
    Idn(1:nLev-1) = cL(:)             + dL(:)*rL(:)*gL(:)          + Y(:)*ExMu0(1:nLev-1)
    Idn(nLev)     = cL(nLay)*gL(nLay) + dL(nLay)*rL(nLay)          + Y(nLay)*ExMu0(nLev)

    if(debug)then
       open(unit=1,file='Idn.out')
       write(1,*) Idn
       close(1)
       open(unit=2,file='Iup.out')
       write(2,*) Iup
       close(2)
       open(unit=2,file='tauL.out')
       write(2,*) tauL 
       close(2)
    end if

    deallocate(aT )
    deallocate(bT )
    deallocate(cT )
    deallocate(dT )
    deallocate(xT )

    return

  end subroutine get_Iup_Idn

  ! ==================================================================================

  subroutine thomas_solve(n,a,b,c,d,x)

    implicit none

    integer, intent(in)    :: n 
    real(8), intent(in)    :: a(n)
    real(8), intent(inout) :: b(n)
    real(8), intent(in)    :: c(n)
    real(8), intent(inout) :: d(n)
    real(8), intent(out)   :: x(n)

    integer k 
    real(8) m 

    ! a_i*x_{i-1} + b_i*x_i + c_i*x_{i+1} = d_i
    !
    ! O[n] efficient, nice.
    ! by default, b and d arrays get written over in process.

    ! forward substitution
    do k=2,n
       m    = a(k)/b(k-1)
       b(k) = b(k) - m*c(k-1)
       d(k) = d(k) - m*d(k-1)
    end do

    ! backward substitution
    x(n) = d(n)/b(n)

    do k=n-1,1,-1
       x(k) = (d(k) - c(k)*x(k+1))/b(k)
    end do

    return

  end subroutine thomas_solve

  ! ==================================================================================

end module mult_scat_core

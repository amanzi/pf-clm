module grid_type_module

  !=========================================================================
  !
  !  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
  !  L                        M  available land surface process model.  
  !  M --COMMON LAND MODEL--  C  
  !  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
  !  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
  !
  !=========================================================================
  ! DESCRIPTION:
  !  Module for grid space variable specification.
  !
  ! REVISION HISTORY:
  !  15 Jan 2000: Paul Houser; Initial code
  !=========================================================================

  use clm_precision
  use clm_infnan
  use clm1d_varpar, only : nlevsoi
  implicit none

  public :: grid_create, grid_init, grid_destroy

  interface grid_create
     module procedure grid_create_1d
     module procedure grid_create_2d
  end interface grid_create

  type,public :: grid_type
     !=== GRID SPACE User-Defined Parameters ==================================

     ! Leaf constants

     real(r8) :: dewmx     ! Maximum allowed dew [mm]

     ! Roughness lengths

     real(r8) :: zlnd      ! Roughness length for soil [m]
     real(r8) :: zsno      ! Roughness length for snow [m]
     real(r8) :: csoilc    ! Drag coefficient for soil under canopy [-]

     ! Soil parameters

     real(r8) :: wtfact    ! Fraction of model area with high water table
     real(r8) :: trsmx0    ! Max transpiration for moist soil+100% veg. [mm/s]
     real(r8) :: scalez    ! Soil layer thickness discretization (m)
     real(r8) :: hkdepth   ! length scale for Ksat decrease(m)

     ! Land surface parameters

     real(r8) :: latdeg    ! Latitude in Degrees 
     real(r8) :: londeg    ! Longitude in Degrees

     real(r8) :: sand(nlevsoi) ! Percent sand in soil 
     real(r8) :: clay(nlevsoi) ! Percent clay in soil 

     real(r8), pointer :: fgrd(:) !Fraction of vegetation class in grid     
     integer , pointer :: pveg(:) !Predominance of vegetation class in grid

     integer :: mask      ! Land=1, Not Land (i.e. not modeled)=0

     !=== CLM Forcing parameters

     real(r8) :: forc_hgt_u ! Observational height of wind [m]
     real(r8) :: forc_hgt_t ! Observational height of temperature [m]
     real(r8) :: forc_hgt_q ! Observational height of humidity [m] 

     !=== Land Surface Fluxes

     real(r8) :: qflx_evap_tot   ! evapotranspiration from canopy height to atmosphere [mm/s]
     real(r8) :: eflx_sh_tot     ! sensible heat from canopy height to atmosphere [W/m2]
     real(r8) :: eflx_lh_tot     ! latent heat flux from canopy height to atmosphere [W/2]
     real(r8) :: eflx_lwrad_out  ! outgoing long-wave radiation from ground+canopy
     real(r8) :: t_ref2m         ! 2 m height air temperature [K]
     real(r8) :: t_rad           ! radiative temperature [K]

     !=== CLM Vegetation parameters

     real(r8) :: rootfr            ! Root Fraction (depth average)

     !=== CLM Soil parameters

     real(r8) :: smpmax       ! Wilting point potential in mm
     integer  :: isoicol      ! Soil color index

     !=== Numerical finite-difference

     real(r8) :: capr         ! Tuning factor to turn first layer T into surface T
     real(r8) :: cnfac        ! Crank Nicholson factor between 0 and 1
     real(r8) :: smpmin       ! Restriction for min of soil poten. (mm)
     real(r8) :: ssi          ! Irreducible water saturation of snow
     real(r8) :: wimp         ! Water impremeable if porosity < wimp
     real(r8) :: pondmx       ! Ponding depth (mm)

     integer  :: tilei        ! Tile index at x,y (or c,r); used to convert tile data into grid data;
     ! works only if there is a single land cover type per grid cell!!

     !=== End Variable List ===================================================
  end type grid_type

contains
  
  function grid_create_2d(ncol, nrow, nt) result(grid)
    type(grid_type),pointer,dimension(:,:) :: grid
    integer, intent(in) :: nrow, ncol, nt

    ! locals
    integer :: r,c

    allocate(grid(1:nrow,1:ncol))
    do c=1,ncol
    do r=1,nrow
       call grid_init(grid(c,r), nt)
    end do
    end do
  end function grid_create_2d

  function grid_create_1d(n, nt) result(grid)
    type(grid_type),pointer,dimension(:,:) :: grid
    integer, intent(in) :: n, nt

    ! locals
    integer :: lcv

    allocate(grid(1:n,1))
    do lcv=1,n
       call grid_init(grid(lcv,1), nt)
    end do
  end function grid_create_1d

  subroutine grid_init(grid, nt)
    type(grid_type) :: grid
    integer, intent(in) :: nt

    grid%dewmx = NaN
    grid%zlnd = NaN
    grid%zsno = NaN
    grid%csoilc = NaN
    grid%wtfact = NaN
    grid%trsmx0 = NaN
    grid%scalez = NaN
    grid%hkdepth = NaN
    grid%latdeg = NaN
    grid%londeg = NaN
    grid%sand(:) = NaN
    grid%clay(:) = NaN

    allocate(grid%fgrd(nt))
    grid%fgrd(:) = NaN
    allocate(grid%pveg(nt))
    grid%pveg(:) = -1

    grid%mask = -1

    grid%forc_hgt_u = NaN
    grid%forc_hgt_t = NaN
    grid%forc_hgt_q = NaN
    grid%qflx_evap_tot = NaN
    grid%eflx_sh_tot = NaN
    grid%eflx_lh_tot = NaN
    grid%eflx_lwrad_out = NaN
    grid%t_ref2m = NaN
    grid%t_rad = NaN

    grid%rootfr = NaN
    grid%smpmax = NaN
    grid%isoicol = -1

    grid%capr = NaN
    grid%cnfac = NaN
    grid%smpmin = NaN
    grid%ssi = NaN
    grid%wimp = NaN
    grid%pondmx = NaN

    grid%tilei = -1
    return
  end subroutine grid_init


  subroutine grid_destroy(grid)
    type(grid_type) :: grid
    if (associated(grid%fgrd)) deallocate(grid%fgrd)
    if (associated(grid%pveg)) deallocate(grid%pveg)
    return
  end subroutine grid_destroy

end module grid_type_module





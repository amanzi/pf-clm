!#include <misc.h>

module tile_type_module
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
!  Module for tile space variable specification.
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!=========================================================================
! $Id: drv_tilemodule.F90,v 1.1.1.1 2006/02/14 23:05:52 kollet Exp $
!=========================================================================

  use clm_precision
  use clm_infnan
  use clm1d_varpar, only : nlevsoi
  implicit none

  public :: tile_create, tile_create_n, tile_init, tile_destroy
  
  type, public:: tile_type
!=== TILE SPACE User-Defined Parameters ====================================
     integer  :: col           ! Grid Column of Tile
     integer  :: row           ! Grid Row of Tile
     integer  :: vegt          ! Vegetation Type of Tile
     integer  :: pveg          ! Predominance of vegetation clas in grid
     real(r8) :: fgrd          ! Fraction of grid covered by a given veg type (%/100)

     real(r8) :: sand(nlevsoi) ! Percent sand in soil (vertically average)
     real(r8) :: clay(nlevsoi) ! Percent clay in soil (vertically average)

     real(r8) :: scalez        ! Soil layer thickness discretization (m)
     real(r8) :: hkdepth       ! length scale for Ksat decrease(m)
     real(r8) :: roota         ! Temporary CLM vegetation parameter
     real(r8) :: rootb         ! Temporary CLM vegetation parameter

!=== End Variable List ===================================================
  end type tile_type

contains
  
  function tile_create() result(tile)
    type(tile_type),pointer :: tile
    allocate(tile)
    call tile_init(tile)
    return
  end function tile_create

  function tile_create_n(n) result(tiles)
    type(tile_type),pointer,dimension(:) :: tiles
    integer, intent(in) :: n

    ! locals
    integer lcv
    
    allocate(tiles(1:n))
    do lcv=1,n
       call tile_init(tiles(lcv))
    end do
    return
  end function tile_create_n

  subroutine tile_init(tile)
    type(tile_type) :: tile
    tile%col = -1
    tile%row = -1
    tile%vegt = -1
    tile%pveg = -1
    tile%fgrd = NaN

    tile%sand(:) = NaN
    tile%clay(:) = NaN

    tile%scalez = NaN
    tile%hkdepth = NaN
    tile%roota = NaN
    tile%rootb = NaN
    return
  end subroutine tile_init

  subroutine tile_destroy(tile)
    type(tile_type) :: tile
    return
  end subroutine tile_destroy
  
end module tile_type_module

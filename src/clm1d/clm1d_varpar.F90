!#include <misc.h>

module clm1d_varpar

!----------------------------------------------------------------------- 
! 
! Purpose: 
! land surface model parameters
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

  use clm_precision
  implicit none

! define level parameters

  integer, parameter :: nlevsoi     =  100   !number of soil levels
  integer, parameter :: nlevlak     =  100   !number of lake levels
  integer, parameter :: nlevsno     =  5    !number of maximum snow levels

  integer, parameter :: numrad      =   2   !number of solar radiation bands: vis, nir
  integer, parameter :: numcol      =   8   !number of soil color types

  integer, parameter :: numtypes    = 18    !number of soil types
end module clm1d_varpar

module clm_precision
!-------------------------------------------------------------------------------
! Purpose:
! Define the precision to use for floating point operations
!-------------------------------------------------------------------------------
  implicit none

!#if (defined DOUBLE_PRECISION)
!  integer, parameter :: r8 = selected_real_kind(12)
!#else
  integer, parameter :: r8 = selected_real_kind(8)
!#endif

  integer,parameter :: MAX_FILENAME_LENGTH = 100
  
end module clm_precision

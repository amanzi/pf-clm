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

end module clm_precision

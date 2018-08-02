module clm_precision
!-------------------------------------------------------------------------------
! Purpose:
! Define the precision to use for floating point operations
!-------------------------------------------------------------------------------
  use ISO_C_BINDING, only: C_DOUBLE,C_INT
  implicit none

!#if (defined DOUBLE_PRECISION)
!  integer, parameter :: r8 = selected_real_kind(12)
!#else

!  integer, parameter :: r8 = selected_real_kind(8)
  integer, parameter :: r8 = C_DOUBLE

!  integer, parameter :: i4 = selected_int_kind(8)
  integer, parameter :: i4 = C_INT
!#endif

end module clm_precision

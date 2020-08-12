module clm_infnan
!-------------------------------------------------------------------------
! Purpose:
!	Set parameters for the floating point flags "inf" Infinity
!	and "nan" not-a-number. As well as "bigint" the point
!	at which integers start to overflow. These values are used
!	to initialize arrays with as a way to detect if arrays
!	are being used before being set.
!-------------------------------------------------------------------------
  use clm_precision
  implicit none
!#if (defined DOUBLE_PRECISION)
 ! real(r8), parameter :: inf = O'777600000000000000000' 
  !real(r8), parameter :: nan = O'777677777777777777777' 
!#else 
  integer :: int_infinity
  integer :: int_nan
  data int_infinity/O'17740000000'/
  data int_nan/O'17757777777'/
  real(r8) inf
  real(r8) nan
  EQUIVALENCE (inf,int_infinity)
  EQUIVALENCE (nan,int_nan)

!  real(r8), parameter :: inf = O'17740000000'
!  real(r8), parameter :: nan = O'17757777777'
  real(r8), parameter :: invalid = -9999.0
!#endif

  integer, parameter  :: bigint = 100000000
end module clm_infnan



  




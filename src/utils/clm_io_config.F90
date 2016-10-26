!default configuration, change me at will
module clm_io_config
  use clm_precision, only : MAX_FILENAME_LENGTH
  implicit none

  character(LEN=50) :: ACCESS='stream'
  character(LEN=50) :: FORM='unformatted'

  public :: io_get_stream

contains

  function io_get_stream() result(iunit)
    integer :: iunit
    logical :: op
    integer,parameter :: maxunit=100
    integer,parameter :: minunit=10

    do iunit=minunit,maxunit
       inquire(unit=iunit, opened=op)
       if (.not. op) then
          return
       endif
    end do

    print *, 'ERROR: all file unit numbers used'
    stop
  end function io_get_stream
  
end module clm_io_config

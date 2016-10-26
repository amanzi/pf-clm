module io_type_module 
  implicit none
  private
  
  ! POD container for io control
  type, public :: io_type
     character(LEN=MAX_FILENAME_LENGTH) :: log_filename         ! name of the logfile
     character(LEN=MAX_FILENAME_LENGTH) :: ranked_log_filename         ! name of the logfile
     character(LEN=MAX_FILENAME_LENGTH) :: output_dir

     integer :: log
     integer :: ranked_log
     integer :: restart_last
     integer :: restart_daily
     integer :: output_in_subdirectories
  end type io_type

  public :: io_open, io_close
  
contains

  
end module io_type_module

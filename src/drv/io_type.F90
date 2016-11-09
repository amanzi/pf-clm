module io_type_module 
  use clm_io_config, only : MAX_FILENAME_LENGTH, io_get_stream
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
     integer :: output_1d
     integer :: dump_interval
     integer :: dump_current
     integer :: output_in_subdirectories
     integer :: write_binaries
  end type io_type

  public :: io_init, io_open, io_close
  
contains

  subroutine io_init(io)
    implicit none
    type(io_type) :: io

    io%log = 0
    io%ranked_log = 0
    io%restart_last = 1
    io%restart_daily = 0
    io%output_1d = 0
    io%dump_interval = -1
    io%dump_current = 0
    io%output_in_subdirectories = 0
    io%write_binaries = 0
  end subroutine io_init
  
  subroutine io_open(io, output_dir, rank, log_by_rank)
    implicit none
    type(io_type) :: io
    character(LEN=MAX_FILENAME_LENGTH),intent(in) :: output_dir
    integer,intent(in) :: rank
    integer,intent(in) :: log_by_rank
    character(len=100) :: rank_string

    if (len_trim(output_dir) > 0) then
       write(io%output_dir,*) output_dir
    end if

    
    if (log_by_rank /= 0) then
       io%ranked_log = io_get_stream()
       open(io%ranked_log, file="clm."//trim(adjustl(rank_string))//".log", action="write")
    else
       io%ranked_log = 0
    endif

    if (rank == 0) then
       io%log = io_get_stream()
       open(io%log, file="clm.log", action="write")
    else
       io%log = 0
    end if
  end subroutine io_open

  subroutine io_close(io)
    implicit none
    type(io_type) :: io

    if (io%ranked_log /= 0) then
       close(io%ranked_log)
    end if

    if (io%log /= 0) then
       close(io%log)
    end if    
  end subroutine io_close
  
  
end module io_type_module

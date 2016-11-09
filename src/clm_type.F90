module clm_type_module
  use clm_precision, only : r8
  use drv_type_module
  use io_type_module
  use tile_type_module
  use grid_type_module
  use clm1d_type_module
  use clm_host
  implicit none
  
  private

  public :: clm_init, &
       clm_setup_begin, &
       clm_setup_end, &
       clm_restart, &
       clm_advance_time, &
       clm_destroy

  !
  ! basic class declaration for all clm data
  ! ----------------------------------------------------
  type, public:: clm_type
     type(drv_type) :: drv
     type(io_type) :: io

     type(tile_type), pointer, dimension(:) :: tile
     type(clm1d_type), pointer, dimension(:) :: clm
     integer :: ntiles
     
     type(grid_type), pointer, dimension(:,:) :: grid
     integer :: grid_nrows, grid_ncols

     integer :: rank

     integer :: istep
     real(r8) :: time
     real(r8) :: dt
  end type clm_type

contains

  
  !
  ! initializer, touches all memory
  !------------------------------------------------------
  subroutine clm_init(this, rank, ncols, nrows, ntypes)
    implicit none
    type(clm_type) :: this
    integer,intent(in) :: rank
    integer,intent(in) :: nrows,ncols,ntypes

    ! set this data
    this%ntiles = nrows * ncols
    this%grid_nrows = nrows
    this%grid_ncols = ncols
    this%rank = rank

    this%istep = 0
    this%dt = 0.
    this%time = 0.

    ! also set some drv data which replicates this data
    call drv_init(this%drv)
    this%drv%nr = this%grid_nrows
    this%drv%nc = this%grid_ncols    

    ! initialize io info
    call io_init(this%io)

    ! nullify for now, need extra arguments to create
    nullify(this%clm)
    this%tile => tile_create_n(this%ntiles) ! note this may be too many
    write(*,*) "Creating n tiles : ", this%ntiles
    this%grid => grid_create_2d(ncols, nrows, ntypes)
    return
  end subroutine clm_init

  !
  ! setup, moves data around assuming all have been init'd
  !------------------------------------------------------
  subroutine clm_setup_begin(this)
    implicit none
    type(clm_type) :: this

    ! locals
    integer :: t

    call drv_g2tile(this%drv, this%grid, this%tile, this%clm, this%ntiles)
    this%ntiles = this%drv%nch ! note this call calculates the actual number of tiles
  end subroutine clm_setup_begin


  !
  ! setup, moves data around assuming all have been init'd
  !------------------------------------------------------
  subroutine clm_setup_end(this)
    implicit none
    type(clm_type) :: this

    ! locals
    integer :: t

    call drv_readvegpf(this%drv, this%grid, this%tile, this%clm)
    do t=1,this%ntiles
       this%clm(t)%istep = this%istep
       call drv_g2clm(this%drv, this%grid, this%tile(t), this%clm(t))
       call drv_clmini(this%drv, this%grid, this%tile(t), this%clm(t))
    end do
  end subroutine clm_setup_end
  
  !
  ! setter for time info
  !------------------------------------------------------------------
  subroutine clm_advance_time(this, host, istep, time, dt)
    implicit none
    type(clm_type) :: this
    type(host_type),intent(in) :: host
    integer,intent(in) :: istep
    real(r8),intent(in) :: dt, time

    ! locals
    integer :: t
    integer :: log_1dout
    
    ! set local time info for logging
    this%istep = istep
    this%time = time
    this%dt = dt

    ! set column time info for doing the work
    this%clm(:)%dtime = dt

    ! tick the driver (sets minutes, seconds, etc)
    if (dt < 1) then ! FIXME --etc
       write(*,*) "CLM cannot take timesteps of smaller than a second"
       stop
    end if
       
    this%drv%ts = nint(dt)
    this%drv%endtime = 0
    call drv_tick(this%drv)

    ! advance the physics
    do t=1,this%ntiles
       if (host%planar_mask(t) == 1) then
          call clm1d_main(this%clm(t), this%drv%day, this%drv%gmt)
       end if
    end do

    ! write 1d output
    if (this%io%output_1d == 1) then
       if (this%io%ranked_log == 0) then
          log_1dout = 0
       else
          log_1dout = 1
       end if
       call drv_1dout(this%drv, this%tile, this%clm, log_1dout)
    end if

  end subroutine clm_advance_time
  

  !
  ! read or write a restart file
  !
  !  1 = read
  !  2 = write
  !------------------------------------------------------
  subroutine clm_restart(rw, this, istep)
    implicit none
    type(clm_type) :: this
    integer, intent(in) :: rw
    integer, intent(in) :: istep

    if (istep == -1) then
       call drv_restart(rw, this%drv, this%tile, this%clm, this%rank, this%istep)
    else
       call drv_restart(rw, this%drv, this%tile, this%clm, this%rank, istep)
    end if
  end subroutine clm_restart
  
  
  !
  ! destructor, cleans memory
  !------------------------------------------------------
  subroutine clm_destroy(this)
    implicit none
    type(clm_type) this

    ! locals
    integer :: t, r, c
    
    call drv_destroy(this%drv)

    if (associated(this%tile)) then
       do t=1,this%ntiles
          call tile_destroy(this%tile(t))
       end do
       deallocate(this%tile)
    end if
          
    if (associated(this%grid)) then
       do r=1,this%grid_nrows
       do c=1,this%grid_ncols
          call grid_destroy(this%grid(r,c))
       end do
       end do
       deallocate(this%grid)
    end if
     
    if (associated(this%clm)) then
       do t=1,this%ntiles
          call clm1d_destroy(this%clm(t))
       end do
       deallocate(this%clm)
    end if
  end subroutine clm_destroy

  


end module clm_type_module

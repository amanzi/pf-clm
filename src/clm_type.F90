module clm_type_module
  use drv_type_module, only : drv_type
  use io_type_module, only : io_type
  use tile_type_module, only : tile_type
  use grid_type_module, only : grid_type
  use clm1d_type_module, only : clm1d_type
  implicit none
  
  private

  public :: clm_create, &
       clm_init, &
       clm_destroy

  !
  ! basic class declaration for all clm data
  ! ----------------------------------------------------
  type, public:: clm_type
     ! type for SAVE data from original subroutine
     type (drv_type) :: drv
     type (io_type) :: io

     type (tile_type), pointer, dimension(:) :: tile
     type (clm1d_type), pointer, dimension(:) :: clm
     integer ntiles
     
     type (grid_type), pointer, dimension(:,:) :: grid
     integer grid_nrows, grid_ncols

     integer rank
  end type clm_type

contains

  !
  ! creates and allocates data structures
  !------------------------------------------------------
  function clm_create(rank) result(this)
    type(clm_type),pointer:: this
    integer :: rank

    allocate(this)
    call clm_init(this, rank)
    return
  end function clm_create

  
  !
  ! initializer, touches all memory
  !------------------------------------------------------
  subroutine clm_init(this, rank)
    type(clm_type) :: this
    integer :: rank

    this%ntiles = -1
    this%grid_nrows = -1
    this%grid_ncols = -1
    
    this%rank = rank
    call drv_init(this%drv)
    call io_init(this%io, rank)

    nullify(this%tile)
    nullify(this%clm)
    nullify(this%grid)

    return
  end subroutine clm_init

  
  !
  ! destructor, cleans memory
  !------------------------------------------------------
  subroutine clm_destroy(this)
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

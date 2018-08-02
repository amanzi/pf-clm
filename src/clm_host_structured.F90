! ------------------------------------------------------------------------------
! CLM Host file
!
! This file, written by the user, provides information on how to map
! data from the host's data structures to CLM's data structures.  Once
! this file is written, following this format, clm_host_transfer
! functions use this capability to do the actual mapping.
!
! Required data type:
!
!  type host_type
!    integer :: ncells          ! number of 3D soil cells
!    integer :: ncells_g        ! length of 3D cell-based arrays
!    integer :: ncolumns        ! number of 2D (map-view) columns
!    integer :: ncolumns_g      ! length of 2D surface-based arrays
!    integer :: planar_mask(1:ncolumns)
!                               ! mask of columns, 1 if column is used, 0 if not
!  end type
!
! Required methods:
!
!  host_cell_index
!  ----------------
!  Maps a CLM cell into any cell-based array from the host code.
!
!  l = host_cell_index(host, i, j, k)
!     integer :: i      ! column index of the clm grid
!     integer :: j      ! row index of the clm grid
!     integer :: k      ! z index of the column of cells
!     type(host_type)   ! host
!
!     integer :: l      ! index into the host code's cell-based array
!
!
!  host_column_index
!  ------------------
!  Maps a CLM cell into any cell-based array from the host code.
!
!  l = host_cell_index(host, i, j)
!     integer :: i      ! column index of the clm grid
!     integer :: j      ! row index of the clm grid
!     type(host_type)   ! host
!
!     integer :: l      ! index into the host code's column-based array
!
!
!
! ParFlow host code (a structured example)
! ============================================
! ParFlow is a logically-structured mesh which uses masks to hide columns.
!
! Indexing is provided through strided vector access.  Most structured
! codes will implement something very similar to this.
!
! Author: Ethan Coon (ecoon _at_ lanl.gov)
! ------------------------------------------------------------------------------

module clm_host
  use clm_precision
  implicit none

  private

  type, public :: host_type
     ! required data
     integer :: ncells
     integer :: ncells_g
     integer :: ncolumns
     integer :: ncolumns_g

     integer, allocatable :: planar_mask(:)  ! mask for inclusion/noninclusion of a column

     ! additional data
     integer :: nx, ny, nz                   ! grid size, locally
     integer :: nx_f, ny_f, nz_f             ! grid size with ghosts
     integer :: j_incr, k_incr               ! useful offset sizes
     integer, allocatable :: topo_mask(:,:)  ! mask for topography: m(1) = top of column,
                                             !     m(2) = bottom of clm column,
                                             !     m(3) = bottom of host column
                                             ! this is only used for the below index methods
  end type host_type

  public :: host_init, host_destroy, &
       host_cell_index, host_column_index, &
       host_write_to_log

contains

  !
  ! Set up the host info
  !
  ! ------------------------------------------------------------------
  subroutine host_init(host, nx, ny, nz, mask)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type) :: host
    integer, intent(in) :: nx, ny, nz
    real(r8), intent(in) :: mask((nx+2)*(ny+2)*(nz+2)) ! mask from host code:
                           ! 0:inactive, 1:active

    ! locals
    integer :: i,j,k,l
    integer :: counter, col_id
    
    host%nx = nx
    host%ny = ny
    host%nz = nz

    host%nx_f = nx + 2
    host%ny_f = ny + 2
    host%nz_f = nz + 2

    host%j_incr = host%nx_f
    host%k_incr = host%nx_f * host%ny_f

    host%ncells = nx * ny * nz
    host%ncells_g = host%nx_f * host%ny_f * host%nz_f
    host%ncolumns = nx * ny
    host%ncolumns_g = host%nx_f * host%ny_f * 3 ! NOTE: ParFlow needs ghost cells in z on surface?

    if (.not.allocated(host%topo_mask)) allocate(host%topo_mask(host%ncolumns,3))
    if (.not.allocated(host%planar_mask)) allocate(host%planar_mask(host%ncolumns))         

    !=== Initialize the mask
    !    This is two components: 
    !    1) a x-y mask of 0/1 for inactive/active and 
    !    2) a z/k mask that takes three values 
    !      (1)= top of LS/PF domain 
    !      (2)= top-nlevsoi and 
    !      (3)= the bottom of the LS/PF domain.
    host%topo_mask(:,:) = 0
    host%planar_mask(:) = 0

    col_id = 0
    do j=1,host%ny
       do i=1,host%nx
          counter = 0
          col_id = col_id + 1

          host%topo_mask(col_id,3) = 1

          do k = host%nz, 1, -1 ! loop over z
             l = 1+i + host%nx_f*j + host%nx_f*host%ny_f*k
             if (mask(l) > 0) then
                counter = counter + 1
                if (counter == 1) then 
                   host%topo_mask(col_id,1) = k
                   host%planar_mask(col_id) = 1
                end if
             endif
             if (mask(l) == 0 .and. mask(l+host%k_incr) > 0) host%topo_mask(col_id,3) = k+1
          enddo ! k
          host%topo_mask(col_id, 2) = host%topo_mask(col_id, 1)-nlevsoi
       enddo
    end do
  end subroutine host_init


  !
  ! Destructor
  !
  ! ------------------------------------------------------------------
  subroutine host_destroy(host)
    implicit none
    type(host_type) :: host

    if (allocated(host%topo_mask)) deallocate(host%topo_mask)
    if (allocated(host%planar_mask)) deallocate(host%planar_mask)
  end subroutine host_destroy

  !
  ! Maps i,j,k, where k is in range 1,nlevsoi, into a structured 3D
  ! array with ghost cells.  Note that this respects the host code's
  ! mask for dead, above-topography cells.
  !
  ! ------------------------------------------------------------------
  function host_cell_index(host, i, j, k) result(l)
    implicit none
    type(host_type), intent(in) :: host
    integer, intent(in) :: i,j,k
    integer :: l

    ! local
    integer :: col_id
    col_id = i + (j-1) * host%ny    
    l = 1 + i + host%j_incr*j + host%k_incr*(host%topo_mask(col_id,1)-(k-1))
  end function host_cell_index


  !
  ! Maps i,j into a structured 2D array with ghost cells
  !
  ! ------------------------------------------------------------------
  function host_column_index(host, i, j) result(l)
    implicit none
    type(host_type), intent(in) :: host
    integer, intent(in) :: i,j
    integer :: l

    l = 1+i + j*host%j_incr + host%k_incr ! NOTE: ParFlow needs ghost cells in z on surface?
  end function host_column_index


  !
  ! write host info to log file
  !
  ! ------------------------------------------------------------------
  subroutine host_write_to_log(host, iounit)
    implicit none
    type(host_type), intent(in) :: host
    integer,intent(in) :: iounit

    write(iounit,*) "Host code: Structured"
    write(iounit,*) "local dimensions:"
    write(iounit,*) '  local NX:',host%nx,' NX with ghost:',host%nx_f
    write(iounit,*) '  local NY:',host%ny,' NY with ghost:',host%ny_f
    write(iounit,*) '  local NZ:',host%nz, 'NZ with ghost:',host%nz_f
    write(iounit,*) "total sizes:"
    write(iounit,*) "  ncells:", host%ncells, "with ghost:", host%ncells_g
    write(iounit,*) "  ncolumns:", host%ncolumns, "with ghost:", host%ncolumns_g
  end subroutine host_write_to_log

end module clm_host

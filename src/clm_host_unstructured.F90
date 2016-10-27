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

     integer, allocatable :: topo_mask(:,:)  ! mask for topography: m(1) = top of column,
                                             !     m(2) = bottom of clm column,
                                             !     m(3) = bottom of host column
     integer, allocatable :: planar_mask(:)  ! mask for inclusion/noninclusion of a column
  end type host_type

  public :: host_create, host_destroy, &
       host_cell_index, host_column_index

contains

  !
  ! Set up the host info
  !
  ! col_inds: ncols x 2 array, consisting of pairs of (bottom, top)
  !           cell number for each column
  !------------------------------------------------------------------------------
  subroutine host_create(info, ncells, ncells_g, ncolumns, ncolumns_g, col_inds)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type) :: info
    integer, intent(in) :: ncells, ncells_g, ncolumns, ncolumns_g
    integer, intent(in) :: col_inds(ncolumns, 2)

    ! local
    integer :: col_id
    
    info%ncells = ncells
    info%ncells_g = ncells_g
    info%ncolumns = ncolumns
    info%ncolumns_g = ncolumns_g

    if (.not.allocated(info%topo_mask)) allocate(info%topo_mask(ncolumns,3))
    if (.not.allocated(info%planar_mask)) allocate(info%planar_mask(ncolumns))         

    !=== Initialize the mask
    !    This is two components: 
    !    1) a x-y mask of 0/1 for inactive/active and 
    !    2) a z/k mask that takes three values 
    !      (1)= top of LS/PF domain 
    !      (2)= top-nlevsoi and 
    !      (3)= the bottom of the LS/PF domain.
    info%planar_mask(:) = 1

    do col_id=1,ncolumns
       info%topo_mask(col_id, 1) = col_inds(col_id, 2)
       info%topo_mask(col_id, 2) = col_inds(col_id, 2) - (nlevsoi - 1)
       info%topo_mask(col_id, 3) = col_inds(col_id, 1)
       call ASSERT(info%topo_mask(col_id, 2) <= info%topo_mask(col_id, 3), "Not enough cells in a column, need at least 10")
    end do

  end subroutine host_create


  !
  ! Destructor
  !
  ! ------------------------------------------------------------------
  subroutine host_destroy(info)
    implicit none
    type(host_type) :: info

    if (allocated(info%topo_mask)) deallocate(info%topo_mask)
    if (allocated(info%planar_mask)) deallocate(info%planar_mask)
  end subroutine host_destroy

  !
  ! Maps i,j,k, where k is in range 1,nlevsoi, into a structured 3D
  ! array with ghost cells.  Note that this respects the host code's
  ! mask for dead, above-topography cells.
  !
  ! ------------------------------------------------------------------
  function host_cell_index(info, i, j, k) result(l)
    implicit none
    type(host_type), intent(in) :: info
    integer, intent(in) :: i,j,k
    integer :: l

    l = info%topo_mask(i,1) - (k-1)
  end function host_cell_index


  !
  ! Maps i,j into a structured 2D array with ghost cells
  !
  ! ------------------------------------------------------------------
  function host_column_index(info, i, j) result(l)
    implicit none
    type(host_type), intent(in) :: info
    integer, intent(in) :: i,j
    integer :: l

    l = i
  end function host_column_index

end module clm_host

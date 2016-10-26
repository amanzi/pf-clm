module clm_host_info
  use clm_precision
  implicit none

  private

  type, public :: host_info_type
     ! required data
     integer :: ncells
     integer :: ncells_g
     integer :: ncolumns
     integer :: ncolumns_g

     integer, allocatable :: planar_mask(:)  ! mask for inclusion/noninclusion of a column
     real(r8), allocatable :: dz(:)          ! grid spacing of each cell

     ! additional data
     integer :: nx, ny, nz                   ! grid size, locally
     integer :: nx_f, ny_f, nz_f             ! grid size with ghosts
     integer :: j_incr, k_incr               ! useful offset sizes
     integer, allocatable :: topo_mask(:,:)  ! mask for topography: m(1) = top of column,
                                             !     m(2) = bottom of clm column,
                                             !     m(3) = bottom of host column
                                             ! this is only used for the below index methods
  end type host_info_type

  public :: host_info_create, host_info_destroy, &
       host_info_cell_index, host_info_column_index

contains

  !
  ! Set up the host info
  !
  ! ------------------------------------------------------------------
  subroutine host_info_create(info, nx, ny, nz, dz, mask)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_info_type) :: info
    integer, intent(in) :: nx, ny, nz
    real(r8), dimension(:), intent(in) :: dz
    real(r8), intent(in) :: mask((nx+2)*(ny+2)*(nz+2)) ! mask from host code:
                           ! 0:inactive, 1:active

    ! locals
    integer :: i,j,k,l
    integer :: counter, col_id
    
    info%nx = nx
    info%ny = ny
    info%nz = nz

    info%nx_f = nx + 2
    info%ny_f = ny + 2
    info%nz_f = nz + 2

    info%j_incr = info%nx_f
    info%k_incr = info%nx_f * info%ny_f

    info%ncells = nx * ny * nz
    info%ncells_g = info%nx_f * info%ny_f * info%nz_f
    info%ncolumns = nx * ny
    info%ncolumns_g = info%nx_f * info%ny_f * 3 ! NOTE: ParFlow needs ghost cells in z on surface?

    if (.not.allocated(info%dz)) allocate(info%dz(info%ncells_g))
    if (.not.allocated(info%topo_mask)) allocate(info%topo_mask(info%ncolumns,3))
    if (.not.allocated(info%planar_mask)) allocate(info%planar_mask(info%ncolumns))         

    !=== Initialize the mask
    !    This is two components: 
    !    1) a x-y mask of 0/1 for inactive/active and 
    !    2) a z/k mask that takes three values 
    !      (1)= top of LS/PF domain 
    !      (2)= top-nlevsoi and 
    !      (3)= the bottom of the LS/PF domain.
    col_id = 0
    do j=1,info%ny
       do i=1,info%nx
          counter = 0
          col_id = col_id + 1

          info%topo_mask(col_id,3) = 1

          do k = info%nz, 1, -1 ! loop over z
             l = 1+i + info%nx_f*j + info%nx_f*info%ny_f*k
             if (mask(l) > 0) then
                counter = counter + 1
                if (counter == 1) then 
                   info%topo_mask(col_id,1) = k
                   info%planar_mask(col_id) = 1
                end if
             endif
             if (mask(l) == 0 .and. mask(l+info%k_incr) > 0) info%topo_mask(col_id,3) = k+1
          enddo ! k
          info%topo_mask(col_id, 2) = info%topo_mask(col_id, 1)-nlevsoi
       enddo
    end do

    !=== Initialize the grid z centroid, dz
    info%dz(:) = dz(:)    
  end subroutine host_info_create


  !
  ! Destructor
  !
  ! ------------------------------------------------------------------
  subroutine host_info_destroy(info)
    implicit none
    type(host_info_type) :: info

    if (allocated(info%dz)) deallocate(info%dz)
    if (allocated(info%topo_mask)) deallocate(info%topo_mask)
    if (allocated(info%planar_mask)) deallocate(info%planar_mask)
  end subroutine host_info_destroy

  !
  ! Maps i,j,k, where k is in range 1,nlevsoi, into a structured 3D
  ! array with ghost cells.  Note that this respects the host code's
  ! mask for dead, above-topography cells.
  !
  ! ------------------------------------------------------------------
  function host_info_cell_index(info, i, j, k) result(l)
    implicit none
    type(host_info_type), intent(in) :: info
    integer, intent(in) :: i,j,k
    integer :: l

    ! local
    integer :: col_id
    col_id = i + j * info%ny    
    l = 1+i + info%j_incr*j + info%k_incr*(info%topo_mask(col_id,1)-(k-1))
  end function host_info_cell_index


  !
  ! Maps i,j into a structured 2D array with ghost cells
  !
  ! ------------------------------------------------------------------
  function host_info_column_index(info, i, j) result(l)
    implicit none
    type(host_info_type), intent(in) :: info
    integer, intent(in) :: i,j
    integer :: l

    l = 1+i + j*info%j_incr + info%k_incr ! NOTE: ParFlow needs ghost cells in z on surface?
  end function host_info_column_index
  
  

  

end module clm_host_info

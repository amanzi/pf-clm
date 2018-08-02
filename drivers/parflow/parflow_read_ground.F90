subroutine parflow_read_ground(host,drv,ix,iy,gnx,gny,latlon,sand,clay,color_index,fractional_ground )

  !=========================================================================
  !
  !  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
  !  L                        M  available land surface process model.  
  !  M --COMMON LAND MODEL--  C  
  !  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
  !  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
  !
  !=========================================================================
  ! DESCRIPTION:
  !  This primary goal of this routine is to determine tile space.
  !
  ! REVISION HISTORY:
  !  15 Jan 2000: Paul Houser; Initial code
  !=========================================================================      
  ! $Id: drv_readvegtf.F90,v 1.1.1.1 2006/02/14 23:05:52 kollet Exp $
  !=========================================================================

  use clm_precision, only : r8
  use drv_type_module, only : drv_type          ! 1-D Land Model Driver variables
  use clm1d_varpar, only : nlevsoi
  use clm_host
  implicit none

  !=== Arguments ===========================================================

  type(host_type),intent(in) :: host
  type(drv_type),intent(in)  :: drv              

  integer,intent(in) :: ix,iy,gnx, gny   !global grid indicies from ParFlow

  real(r8),intent(out) :: latlon(2, host%ncolumns_g)    ! latitude,longitude [degrees]
  real(r8),intent(out) :: sand(host%ncells_g)          ! percent sand FIXME: 0-1 or 0-100? --etc
  real(r8),intent(out) :: clay(host%ncells_g)          ! percent clay FIXME: 0-1 or 0-100? --etc
  integer,intent(out) :: color_index(host%ncolumns_g)  ! color index FIXME: document! --etc
  real(r8),intent(out) :: fractional_ground(drv%nt, host%ncolumns_g) ! fraction of land surface of type t
  
  !=== Local Variables =====================================================

  integer  :: c,r,i,j,k,l,m     !Loop counters
  real(r8) :: sand_tmp        !temporary value of input sand
  real(r8) :: clay_tmp        !temporary value of input clay

  !=== End Variable Definition =============================================

  !=== Read in Vegetation Data
  !open(2,file=trim(adjustl(drv%vegtf))//'.'//trim(adjustl(RI)),form='formatted',action='read')
  open(2,file=trim(adjustl(drv%vegtf)),form='formatted',action='read')

  read(2,*)  !skip header
  read(2,*)  !skip header
  ! do r=1,drv%nr     !rows
   ! do c=1,drv%nc  !columns
  do r =1, gny  ! @RMM replaced local row/column with global grid
     do c = 1, gnx
        if (((c > ix).and.(c <= (ix+host%nx))).and.((r > iy).and.(r <= (iy+host%ny)))) then
           l = host_column_index(host,c-ix,r-iy)
           read(2,*) i,j,          &
                latlon(1, l), &
                latlon(2, l), &
                sand_tmp,              &
                clay_tmp,              &
                color_index(l), &
                (fractional_ground(m,l),m=1,drv%nt)

           do k=1,nlevsoi
              l = host_cell_index(host,c-ix,r-iy,k)
              sand(l) = sand_tmp
              clay(l) = clay_tmp
           end do
        else
           read(2,*)
        end if
     enddo ! C 
  enddo ! R 

  close(2)
  return

end subroutine parflow_read_ground

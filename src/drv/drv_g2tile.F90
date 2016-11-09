subroutine drv_g2tile(drv,grid,tile,clm,ntiles)

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

  use clm_precision
  use drv_type_module          ! 1-D Land Model Driver variables
  use tile_type_module      ! Tile-space variables
  use clm1d_type_module             ! 1-D CLM variables
  use grid_type_module      ! Grid-space variables
  implicit none

  !=== Arguments ===========================================================
  integer,intent(in) :: ntiles
  type (drv_type)  :: drv              
  type (tile_type) :: tile(ntiles)
  type (clm1d_type)   :: clm (ntiles)
  type (grid_type) :: grid(drv%nc,drv%nr)   

  !=== Local Variables =====================================================

  integer  :: c,r,t,i,j     !Loop counters
  real(r8) :: rsum          !Temporary vegetation processing variable
  real(r8) :: fvt(drv%nt)   !Temporary vegetation processing variable
  real(r8) :: max           !Temporary vegetation processing variable

  !=== End Variable Definition =============================================
  do r=1,drv%nr     !rows
     do c=1,drv%nc  !columns
        rsum=0.0
        do t=1,drv%nt
           rsum=rsum+grid(c,r)%fgrd(t)
        enddo
        if (rsum >= drv%mina) then
           grid(c,r)%mask=1
        else
           grid(c,r)%mask=0
        endif
     enddo ! C 
  enddo ! R 

  !=== Exclude tiles with MINA (minimum tile grid area),  
  !=== normalize remaining tiles to 100%
  do r=1,drv%nr  !rows
   do c=1,drv%nc  !columns         

        rsum=0.0
        do t=1,drv%nt
           if (grid(c,r)%fgrd(t).lt.drv%mina) grid(c,r)%fgrd(t) = 0.0    ! impose area percent cutoff
           rsum = rsum + grid(c,r)%fgrd(t)
        enddo

        if (rsum.gt.0.0) then ! Renormalize veg fractions within a grid to 1  
           do t=1,drv%nt      ! Renormalize SUMT back to 1.0
              grid(c,r)%fgrd(t) = grid(c,r)%fgrd(t) / rsum
           enddo
        endif

     enddo
  enddo

  !=== Exclude tiles with MAXT (Maximum Tiles per grid), 
  !=== normalize remaining tiles to 100%
  !=== Determine the grid predominance order of the tiles
  !=== PVEG(NT) will contain the predominance order of tiles

  do r=1,drv%nr  !rows
     do c=1,drv%nc  !columns

        do t=1,drv%nt
           fvt(t)=grid(c,r)%fgrd(t)  !fvt= temp fgrd working array
           grid(c,r)%pveg(t)=0
        enddo
        do i=1,drv%nt  !Loop through predominance level
           max=0.0
           t=0
           do j=1,drv%nt
              if (fvt(j) > max)then
                 if (grid(c,r)%fgrd(j) > 0) then
                    max=fvt(j)
                    t=j
                 endif
              endif
           enddo
           if (t > 0)then
              grid(c,r)%pveg(t)=i
              fvt(t)=-999.0       !eliminate chosen from next search 
           endif
        enddo
     enddo !IR
  enddo !IC 

  !=== Impose MAXT Cutoff

  do r=1,drv%nr  !rows
     do c=1,drv%nc  !columns         
        rsum=0.0
        do t=1,drv%nt
           if (grid(c,r)%pveg(t).lt.1) then
              grid(c,r)%fgrd(t)=0.0    
              grid(c,r)%pveg(t)=0  
           else if (grid(c,r)%pveg(t)>drv%maxt) then
              grid(c,r)%fgrd(t)=0.0              ! impose maxt cutoff
              grid(c,r)%pveg(t)=0  
           endif
           rsum=rsum+grid(c,r)%fgrd(t)
        enddo

        if (rsum > 0.0) then   ! Renormalize veg fractions within a grid to 1
           do t=1,drv%nt       ! Renormalize SUMT back to 1.0
              grid(c,r)%fgrd(t) = grid(c,r)%fgrd(t) / rsum
           enddo
        endif

     enddo
  enddo

  !=== Make Tile Space
  ! FIXME: this is rather bad.  There are at least two ways to mask tiles:
  !  1. through no land (i.e. tile%fgrnd == 0)
  !  2. through planar_mask (i.e. no valid cells)
  ! We allocate nrow * ncol tiles, but then only use ones that satisfy
  ! the above.  Should really only allocate and loop over actual used
  ! tiles. --etc

  drv%nch=0
  do t=1,drv%nt                                              !loop through each tile type
     do r=1,drv%nr                                           !loop through rows
        do c=1,drv%nc                                        !loop through columns
           if (grid(c,r)%mask.eq.1) then                     !we have land 
              if (grid(c,r)%fgrd(t) > 0.0) then
                 drv%nch = drv%nch+1                         !Count the number of tiles
                 grid(c,r)%tilei       = drv%nch             !@ Index to convert tile to grid data in one sweep; works only of 1 l-cover per cell
                 tile(drv%nch)%row     = r                   !keep track of tile row
                 tile(drv%nch)%col     = c                   !keep track of tile column
                 tile(drv%nch)%vegt    = t                   !keep track of tile surface type
                 tile(drv%nch)%fgrd    = grid(c,r)%fgrd(t)   !keep track of tile fraction
                 tile(drv%nch)%pveg    = grid(c,r)%pveg(t)   !Predominance of vegetation class in grid
                 tile(drv%nch)%sand(:) = grid(c,r)%sand(:)   !Percent sand in soil
                 tile(drv%nch)%clay(:) = grid(c,r)%clay(:)   !Percent clay in soil
                 clm(drv%nch)%londeg   = grid(c,r)%londeg    !Longitude of tile (degrees)
                 clm(drv%nch)%latdeg   = grid(c,r)%latdeg    !Latitude of tile (degrees)
                 clm(drv%nch)%isoicol  = grid(c,r)%isoicol   !Soil color 
                 clm(drv%nch)%lat = clm(drv%nch)%latdeg*4.*atan(1.)/180. !tile latitude  (radians)
                 clm(drv%nch)%lon = clm(drv%nch)%londeg*4.*atan(1.)/180. !tile longitude (radians)
              endif
           endif
        enddo
     enddo
  enddo
  return

end subroutine drv_g2tile

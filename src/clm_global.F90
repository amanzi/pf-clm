module clm_global_module
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  use drv_gridmodule      ! Grid-space variables
  use clmtype             ! CLM tile variables
  implicit none

public

  ! type for SAVE data from original subroutine
  type (drvdec)          :: drv
  type (tiledec),pointer :: tile(:)
  type (griddec),pointer :: grid(:,:)
  type (clm1d),pointer   :: clm(:)

end module clm_global_module

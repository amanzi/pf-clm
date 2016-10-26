module clm_global_module
  implicit none

  private

  ! type for SAVE data from original subroutine
  type (drvdec)          :: drv
  type (tiledec),pointer :: tile(:)
  type (griddec),pointer :: grid(:,:)
  type (clm1d),pointer   :: clm(:)

end module clm_global_module

module MLMC
  implicit none

CONTAINS
  ! print statements in this module use # 300-399

  subroutine UQ_MLMC( icase )
  !This subroutine is the main driver for the MLMC method.
  !It's first implementation is planned to be using FEDiffSn as the deterministic solver
  !to demonstrate MLMC.  It's second implementation is planned to be using MC transport
  !as the solver.  Potentially the third is with MC transport over random binary media.
  integer :: icase

  print *,"I'm in UQ_MLMC"

  end subroutine UQ_MLMC


end module MLMC

module MLMC
  implicit none

CONTAINS
  ! print statements in this module use # 300-399

  subroutine UQ_MLMC( icase )
  !This subroutine is the main driver for the MLMC method.
  !It's first implementation is planned to be using FEDiffSn as the deterministic solver
  !to demonstrate MLMC.  It's second implementation is planned to be using MC transport
  !as the solver.  Potentially the third is with MC transport over random binary media.
  use MLMCvars, only: MLMC_TOL, numMLMCcells
  integer :: icase

  integer :: Level
  real(8) :: MLMCerrest
  logical :: flMLMC

  call MLMCallocate                                    !allocate variables
  call MLMCinitialize( flMLMC, Level )                 !initialize variables

  !main MLMC loop
  do while(flMLMC)

    !1 Determine number of cells, make array which holds number of cells as a function of level L
    if(Level /= 0) call MLMCaddLevel( Level )

    !2 Solve M~ and V~ for each level L

    !3 Using V~s, compute optimal M~s

    !4 Evaluate any new samples needed

    !MLMCerrest = calcMLMCerrest()!5 test total error, set flMLMC==.false.?
    MLMCerrest = 0.1d0
    if(Level>1 .and. MLMCerrest<=MLMC_TOL) then
      flMLMC = .false.
    else
      Level = Level + 1
print *,"Level:",Level
if(Level==5) flMLMC = .false.
    endif
  enddo !loops over realizations

  end subroutine UQ_MLMC




  subroutine MLMCallocate
  !Allocate arrays here.  Try to only allocate in order to make paralellization easier later.
  use MLMCvars, only: numMLMCcells

  if(.not.allocated(numMLMCcells)) allocate(numMLMCcells(0:0))

  end subroutine MLMCallocate


  subroutine MLMCinitialize( flMLMC, Level )
  !Initialize/load values to variables here.
  use MLMCvars, only: numMLMCcells, numcellsLevel0
  integer :: Level
  logical :: flMLMC

  integer :: i

  Level  = 0
  flMLMC = .true.  

  numMLMCcells(i) = numcellsLevel0

  end subroutine MLMCinitialize




  subroutine MLMCaddLevel( Level )
  !Increase size of arrays to handle new Level.
  use MLMCvars, only: numMLMCcells, nextLevelFactor
  integer :: Level

  integer, allocatable :: tnumMLMCcells(:)

print *,"numMLMCcells:",numMLMCcells
  !add new Level to numMLMCcells and populate it
  call move_alloc(numMLMCcells,tnumMLMCcells)
  allocate(numMLMCcells(0:size(tnumMLMCcells)))
  numMLMCcells(0:size(tnumMLMCcells)-1) = tnumMLMCcells
  numMLMCcells(Level) = numMLMCcells(Level-1) * nextLevelFactor
print *,"numMLMCcells:",numMLMCcells

  end subroutine MLMCaddLevel



end module MLMC

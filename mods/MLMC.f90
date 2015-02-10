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
    !1 aka add new Level and necessary cells
    if(Level /= 0) call MLMCaddLevel( Level )

    !2 Using M~(current samples) solve V~(estimated variance) for each level L

    !3 Using V~s(estimated variance), compute optimal M~s(estimated opt num of samples)
    !3.2 expand arrays to hold new M~s(number of samples)

    !4 Evaluate any new samples needed

    !5 test total error, set flMLMC==.false.?
    !MLMCerrest = calcMLMCerrest()
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
  use MLMCvars, only: numMLMCcells, numcellsLevel0, MLMC_failprob, C_alpha
  use utilities, only: erfi
  integer :: Level
  logical :: flMLMC

  integer :: i

  Level  = 0
  flMLMC = .true.  

  numMLMCcells(i) = numcellsLevel0

  C_alpha = sqrt(2.0d0)*erfi(1.0d0-MLMC_failprob)

  end subroutine MLMCinitialize




  subroutine MLMCaddLevel( Level )
  !Increase size of arrays to handle new Level.
  use MLMCvars, only: numMLMCcells, nextLevelFactor
  integer :: Level

  integer, allocatable :: tnumMLMCcells(:)

  !add new Level to numMLMCcells and populate it
  call move_alloc(numMLMCcells,tnumMLMCcells)
  allocate(numMLMCcells(0:size(tnumMLMCcells)))
  numMLMCcells(0:size(tnumMLMCcells)-1) = tnumMLMCcells
  numMLMCcells(Level) = numMLMCcells(Level-1) * nextLevelFactor

  end subroutine MLMCaddLevel



end module MLMC

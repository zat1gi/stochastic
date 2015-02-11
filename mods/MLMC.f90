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
    !if(Level /= 0) call MLMCsolveEstVar

    !3 Using V~s(estimated variance), compute optimal M~s(estimated opt num of samples)
    !3.2 expand arrays to hold new M~s(number of samples)
    !if(Level /= 0) call MLMCcomputeOptSamps( Level )
    if(Level /= 0) call MLMCaddSamples( Level )

    !4 Evaluate any new samples needed
    call MLMCevalNewSamps( Level )

    !5 test total error, set flMLMC==.false.?
    !MLMCerrest = MLMCcalcErrEst( Level )
    MLMCerrest = 0.1d0
    if(Level>1 .and. MLMCerrest<=MLMC_TOL) then
      flMLMC = .false.
    else
      Level = Level + 1
print *,"Level:",Level
if(Level==5) flMLMC = .false.
    endif
  enddo !loops over realizations

  call MLMCdeallocate

  end subroutine UQ_MLMC




  subroutine MLMCallocate
  !Allocate arrays here.  Try to only allocate in order to make paralellization easier later.
  use MLMCvars, only: numMLMCcells, M_optsamps, uflux, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, bnumMLMCsamps, numcellsLevel0

  if(.not.allocated(numMLMCcells)) allocate(numMLMCcells(0:0))
  if(.not.allocated(M_optsamps)) allocate(M_optsamps(2,0:0))

  if(.not.allocated(uflux)) allocate(uflux(bnumMLMCsamps,0:0,numcellsLevel0))
  if(.not.allocated(Q_ufunctional)) allocate(Q_ufunctional(bnumMLMCsamps,0:0))
  if(.not.allocated(G_ufunctional)) allocate(G_ufunctional(bnumMLMCsamps,0:0))
  if(.not.allocated(Gave)) allocate(Gave(0:0))
  if(.not.allocated(Gvar)) allocate(Gvar(0:0))

  end subroutine MLMCallocate


  subroutine MLMCinitialize( flMLMC, Level )
  !Initialize/load values to variables here.
  use MLMCvars, only: numMLMCcells, numcellsLevel0, MLMC_failprob, C_alpha, &
                      M_optsamps, uflux, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, bnumMLMCsamps
  use utilities, only: erfi
  integer :: Level
  logical :: flMLMC

  Level  = 0
  flMLMC = .true.  

  C_alpha = sqrt(2.0d0)*erfi(1.0d0-MLMC_failprob)

  numMLMCcells    = numcellsLevel0
  M_optsamps(1,0) = bnumMLMCsamps
  M_optsamps(2,0) = 0

  uflux         = 0.0d0
  Q_ufunctional = 0.0d0
  G_ufunctional = 0.0d0
  Gave          = 0.0d0
  Gvar          = 0.0d0

  end subroutine MLMCinitialize


  subroutine MLMCdeallocate
  !Deallocate arrays here.
  use MLMCvars, only: numMLMCcells, M_optsamps, uflux, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar

  if(allocated(numMLMCcells)) deallocate(numMLMCcells)
  if(allocated(M_optsamps)) deallocate(M_optsamps)

  if(allocated(uflux)) deallocate(uflux)
  if(allocated(Q_ufunctional)) deallocate(Q_ufunctional)
  if(allocated(G_ufunctional)) deallocate(G_ufunctional)
  if(allocated(Gave)) deallocate(Gave)
  if(allocated(Gvar)) deallocate(Gvar)

  end subroutine MLMCdeallocate





  subroutine MLMCaddLevel( Level )
  !Increase size of arrays to handle new Level and in case of response function
  !the associated number of cells as well.
  use MLMCvars, only: numMLMCcells, nextLevelFactor, uflux, Q_ufunctional, &
                      G_ufunctional, Gave, Gvar, M_optsamps, bnumMLMCsamps
  integer :: Level

  integer, allocatable :: tiarray1(:)
  integer, allocatable :: tiarray2(:,:)

  real(8), allocatable :: trarray1(:)
  real(8), allocatable :: trarray2(:,:)
  real(8), allocatable :: trarray3(:,:,:)

  !add new Level to numMLMCcells and populate it
  call move_alloc(numMLMCcells,tiarray1)
  allocate(numMLMCcells(0:size(tiarray1)))
  numMLMCcells = 0
  numMLMCcells(0:size(tiarray1)-1) = tiarray1
  numMLMCcells(Level) = numMLMCcells(Level-1) * nextLevelFactor
  deallocate(tiarray1)

  !add new Level to M_optsamps and populate it
  call move_alloc(M_optsamps,tiarray2)
  allocate(M_optsamps(2,0:size(tiarray2(1,:))))
  M_optsamps = 0
  M_optsamps(:,0:size(tiarray2(1,:))-1) = tiarray2
  M_optsamps(1,Level) = bnumMLMCsamps
  deallocate(tiarray2)


  !add new Level and cells to uflux and initialize them
  call move_alloc(uflux,trarray3)
  allocate(uflux(size(trarray3(:,1,1)),0:size(trarray3(1,:,1)),maxval(numMLMCcells(:))))
  uflux = 0.0d0
  uflux(:,0:size(trarray3(1,:,1))-1,1:size(trarray3(1,1,:))) = trarray3
  deallocate(trarray3)

  !add new Level to Q_ufunctional and initialize it
  call move_alloc(Q_ufunctional,trarray2)
  allocate(Q_ufunctional(size(trarray2(:,1)),0:size(trarray2(1,:))))
  Q_ufunctional = 0.0d0
  Q_ufunctional(:,0:size(trarray2(1,:))-1) = trarray2
  deallocate(trarray2)

  !add new Level to G_ufunctional and initialize it
  call move_alloc(G_ufunctional,trarray2)
  allocate(G_ufunctional(size(trarray2(:,1)),0:size(trarray2(1,:))))
  G_ufunctional = 0.0d0
  G_ufunctional(:,0:size(trarray2(1,:))-1) = trarray2
  deallocate(trarray2)

  !add new Level to Gave and initialize it
  call move_alloc(Gave,trarray1)
  allocate(Gave(0:size(trarray1(:))))
  Gave = 0.0d0
  Gave(0:size(trarray1(:))-1) = trarray1
  deallocate(trarray1)

  !add new Level to Gvar and initialize it
  call move_alloc(Gvar,trarray1)
  allocate(Gvar(0:size(trarray1(:))))
  Gvar = 0.0d0
  Gvar(0:size(trarray1(:))-1) = trarray1
  deallocate(trarray1)

  end subroutine MLMCaddLevel


  subroutine MLMCaddSamples( Level )
  !Increase size of arrays to handle new optimal number of samples.
  use MLMCvars, only: uflux, Q_ufunctional, G_ufunctional, M_optsamps
  integer :: Level

  real(8), allocatable :: trarray2(:,:)
  real(8), allocatable :: trarray3(:,:,:)

  !add space for new samples to uflux and initialize the space
  call move_alloc(uflux,trarray3)
  allocate(uflux(maxval(M_optsamps),0:size(trarray3(1,:,1))-1,size(trarray3(1,1,:))))
  uflux = 0.0d0
  uflux(1:size(trarray3(:,1,1)),:,:) = trarray3
  deallocate(trarray3)

  !add space for new samples to Q_functional and initialize the space
  call move_alloc(Q_ufunctional,trarray2)
  allocate(Q_ufunctional(maxval(M_optsamps),0:size(trarray2(1,:))-1))
  Q_ufunctional = 0.0d0
  Q_ufunctional(1:size(trarray2(:,1)),:) = trarray2
  deallocate(trarray2)

  !add space for new samples to G_ufunctional and initialize the space
  call move_alloc(G_ufunctional,trarray2)
  allocate(G_ufunctional(maxval(M_optsamps),0:size(trarray2(1,:))-1))
  G_ufunctional = 0.0d0
  G_ufunctional(1:size(trarray2(:,1)),:) = trarray2
  deallocate(trarray2)

  end subroutine MLMCaddSamples



  !2 Using M~(current samples) solve V~(estimated variance) for each level L
!  subroutine MLMCsolveEstVar( Level )


!  end subroutine MLMCsolveEstVar


  !4 Evaluate any new samples needed
  subroutine MLMCevalNewSamps( Level )
  !Evaluate any new samples needed, recompute functional values
  use MLMCvars, only: M_optsamps
  use rngvars, only: setrngappnum, rngappnum, rngstride
  use mcnp_random, only: RN_init_particle
  integer :: Level

  integer :: ilevel,isamp, isamplow

  do ilevel = 0,Level-1                                        !search each level
    if( M_optsamps(1,ilevel)>M_optsamps(2,ilevel) ) then       !if new samps to compute
      isamplow = merge(M_optsamps(2,ilevel)+1,1,ilevel/=Level) !set lowest samp num
      do isamp = isamplow,M_optsamps(1,ilevel)                 !cycle through samps to compute
        call setrngappnum('MLMCsamp')                          !set rng unique to sample
        call RN_init_particle( int(rngappnum*rngstride+isamp,8) )
                                                               !update solver input info
                                                               !run solver
                                                               !collect uflux
        !run simulation, interface with FEDiffSn to collect response function and store as part of uflux
      enddo
    endif
    !calc Q_ufunctional
    !calc G_ufunctional
    !calc Gave
    !calc Gvar
  enddo

  end subroutine MLMCevalNewSamps


end module MLMC

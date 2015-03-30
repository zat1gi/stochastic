module MLMC
  implicit none

CONTAINS
  ! print statements in this module use # 300-399

  subroutine UQ_MLMC( icase )
  !This subroutine is the main driver for the MLMC method.
  !It's first implementation is planned to be using FEDiffSn as the deterministic solver
  !to demonstrate MLMC.  It's second implementation is planned to be using MC transport
  !as the solver.  Potentially the third is with MC transport over random binary media.
  use MLMCvars, only: MLMC_TOL, numMLMCcells, def_ufunct, MLMCerrest
  use FEDiffSn, only: setflvarspassedtrue
  integer :: icase
  integer :: ifunct

  integer :: Level
  logical :: flMLMC
  logical :: flread = .false.

  call MLMCallocate                                    !allocate variables
  call MLMCinitialize( flMLMC, Level )                 !initialize variables
  call setflvarspassedtrue                             !tell FE mod to accept input from here

  !main MLMC loop
  do while(flMLMC)
    print *,"Level:",Level

    !1 Determine number of cells, make array which holds number of cells as a function of level L
    !1 aka add new Level and necessary cells
    if(Level /= 0) print *,"--reallocate to accomodate new Level--"
    if(Level /= 0) call MLMCaddLevel( Level )
    if(Level /= 0) print *,"--reallocate to accomodate new Level--"
    if(Level /= 0) print *

    !2 Solve initial samples for Level, solve V~(estimated variance) for each level L
    print *,"--evaluate baseline samples, get variance estimate--"
    call MLMCevalNewSamps( Level,icase )
    print *,"--evaluate baseline samples, get variance estimate--"
    print *

    !3 Using V~s(estimated variance), compute optimal M~s(estimated opt num of samples)
    !3.2 expand arrays to hold new M~s(number of samples)
    print *,"--compute optimal number of samples--"
    call MLMCcomputeOptSamps( Level )
    print *,"--compute optimal number of samples--"
    print *
    print *,"--reallocate to accomodate new samples--"
    call MLMCaddSamples( Level )
    print *,"--reallocate to accomodate new samples--"
    print *

    !4 Evaluate any new samples needed, Gave and Gvar calculated here
    print *,"--evaluate extra samples needed--"
    call MLMCevalNewSamps( Level,icase )
    print *,"--evaluate extra samples needed--"
    print *

    !5 test total error, set flMLMC==.false.?
    if(Level>1) print *,"--calculate estimated error--"
    if(Level>1) then
      do ifunct=1,size(def_ufunct(:,1))                  !calculate estimated error for each functional
        MLMCerrest(ifunct) = MLMCcalcErrEst( Level,ifunct )
        print *,"          MLMCerrest :",MLMCerrest(ifunct)
        if(MLMCerrest(ifunct)>MLMC_TOL) then
          print *,"          MLMC_TOL   :",MLMC_TOL," not converged"
        else
          print *,"          MLMC_TOL   :",MLMC_TOL," converged"
        endif
      enddo
    endif

    do ifunct=1,size(def_ufunct(:,1))                    !find which functional to converge
      if(def_ufunct(ifunct,4)==1) exit
    enddo

    if(Level>1 .and. MLMCerrest(ifunct)<=MLMC_TOL) then  !test for convergence
      flMLMC = .false.
    else
      Level = Level + 1
    endif
    if(Level>1) print *,"--calculate estimated error--"

  enddo !loops over realizations

  call MLMCprintfunctionaldata
  call MLMCdeallocate

  end subroutine UQ_MLMC



  subroutine UQ_spatialconv( icase )
  !This subroutine performs a spatial convergence study for all chosen functionals in an
  !effort to determine the spatial convergence parameter for each problem, solve, and QoI
  use genRealzvars, only: s
  use MLMCvars, only: Q_ufunctional, num_ufunct, spatial_Level, ncellwidth, numcellsLevel0, &
                      nextLevelFactor
  use FEDiffSn, only: setflvarspassedtrue, a
  integer :: icase,ilevel

  if(allocated(Q_ufunctional)) deallocate(Q_ufunctional)
  allocate(Q_ufunctional(num_ufunct,1,0:spatial_Level))
  Q_ufunctional = 0.0d0
  if(allocated(ncellwidth)) deallocate(ncellwidth)
  allocate(ncellwidth(0:spatial_Level))
  ncellwidth = 0.0d0
  call setflvarspassedtrue                             !tell FE mod to accept input from here
  a             = s

  do ilevel=0,spatial_Level
    ncellwidth(ilevel) = s/real(numcellsLevel0*nextLevelFactor**ilevel,8)
    call sampleconvInput( ilevel )          !samples average values
    call solveSamples( ilevel,1,icase )     !solves QoIs, isamp=1
  enddo
  call spatial_calcerr_print                   !calc and print errs for functs
  end subroutine UQ_spatialconv


 subroutine iter_calcerr_print
  !This subroutrine calculates the error (based off the most converged iteration)
  !for each functional and iterative convergence values and prints each to a file.
  use MLMCvars, only: Q_ufunctional, num_ufunct, def_ufunct, &
                      numcellsLevel0, nextLevelFactor

  integer :: iiter, ifunct, ibase, igap, lastiter
  real(8), allocatable :: err_ufunctional(:,:,:),slope

  lastiter = size(Q_ufunctional(1,1,:))-1
  if(.not.allocated(err_ufunctional)) allocate(err_ufunctional(num_ufunct,1,0:lastiter))
  err_ufunctional = 0.0d0

  do iiter=0,lastiter
    do ifunct=1,num_ufunct
      err_ufunctional(ifunct,1,iiter) = &
                  abs( Q_ufunctional(ifunct,1,lastiter)-Q_ufunctional(ifunct,1,iiter) ) / &
                       Q_ufunctional(ifunct,1,lastiter)
    enddo
  enddo
  !Here I calculate the log-log slopes and print in each variation possible
  !for each functional chosen.  I calculate slope between the first and second
  !data points, first and third, etc, then first and third, first and fourth, etc.
  !Each of these slope calculations ought to be about the same, but this is a sanity
  !check that they are!
!  call system("test -e plots/MLMCfuncts/spatialslopes.out && rm plots/MLMCfuncts/spatialslopes.out")

!  open(unit=24, file="plots/MLMCfuncts/spatialslopes.out")

!  1130 format(f15.7)
!  1131 format("slope for ")
!  do ifunct=1,num_ufunct
!    write(24,1131,advance="no")
!    if(def_ufunct(ifunct,3)==1) then
!      write(24,1121) def_ufunct(ifunct,1),def_ufunct(ifunct,2)
!    elseif(def_ufunct(ifunct,3)==2) then
!      write(24,1122) def_ufunct(ifunct,1),def_ufunct(ifunct,2)
!    elseif(def_ufunct(ifunct,3)==3) then
!      write(24,1123) def_ufunct(ifunct,1)
!    endif
!    do igap=1,spatial_Level-1
!      do ibase=0,spatial_Level-igap-1
!        slope = spatiallogslope(err_ufunctional,ifunct,ibase,ibase+igap)
!        write(24,1130,advance="no") slope
!      enddo
!      write(24,1124)
!    enddo
!  enddo

!  close(unit=24)      


  !Here we print the actual functional error data so that it can be examined
  !by hand and/or plotted.
!  call system("test -e plots/MLMCfuncts/spatialconv.out && rm plots/MLMCfuncts/spatialconv.out")

!  open(unit=24, file="plots/MLMCfuncts/spatialconv.out")
!  1120 format("#ilevel     num of cells   ")
!  1121 format("L1 cell",i5," to",i5,"  ")
!  1122 format("L2 cell",i5," to",i5,"  ")
!  1123 format(" center cell",i5,"     ")
!  1124 format(" ")
!  write(24,1120,advance="no")
!  do ifunct=1,num_ufunct
!    if(def_ufunct(ifunct,3)==1) then
!      write(24,1121,advance="no") def_ufunct(ifunct,1),def_ufunct(ifunct,2)
!    elseif(def_ufunct(ifunct,3)==2) then
!      write(24,1122,advance="no") def_ufunct(ifunct,1),def_ufunct(ifunct,2)
!    elseif(def_ufunct(ifunct,3)==3) then
!      write(24,1123,advance="no") def_ufunct(ifunct,1)
!    endif
!  enddo
!  write(24,1124)

!  1125 format(i5,"     ",i13,"    ")
!  1126 format(es14.7,"        ")
!  do ilevel=0,spatial_Level-1
!    write(24,1125,advance="no") ilevel,numcellsLevel0*nextLevelFactor**ilevel
!    do ifunct=1,num_ufunct
!      write(24,1126,advance="no") err_ufunctional(ifunct,1,ilevel)
!    enddo
!    write(24,1124)
!  enddo
!  close(unit=24)    

  end subroutine iter_calcerr_print



  subroutine spatial_calcerr_print
  !This subroutrine calculates the error (based off the most converged level)
  !for each functional and the log-log slopes (convergence rates),
  !and prints each of these values to a file.
  use MLMCvars, only: Q_ufunctional, spatial_Level, num_ufunct, def_ufunct, &
                      numcellsLevel0, nextLevelFactor

  integer :: ilevel, ifunct, ibase, igap
  real(8), allocatable :: err_ufunctional(:,:,:),slope

  if(.not.allocated(err_ufunctional)) allocate(err_ufunctional(num_ufunct,1,0:spatial_Level-1))
  err_ufunctional = 0.0d0

  do ilevel=0,spatial_Level-1
    do ifunct=1,num_ufunct
      err_ufunctional(ifunct,1,ilevel) = &
                  abs( Q_ufunctional(ifunct,1,spatial_Level)-Q_ufunctional(ifunct,1,ilevel) ) / &
                       Q_ufunctional(ifunct,1,spatial_Level)
    enddo
  enddo

  !Here I calculate the log-log slopes and print in each variation possible
  !for each functional chosen.  I calculate slope between the first and second
  !data points, first and third, etc, then first and third, first and fourth, etc.
  !Each of these slope calculations ought to be about the same, but this is a sanity
  !check that they are!
  call system("test -e plots/MLMCfuncts/spatialslopes.out && rm plots/MLMCfuncts/spatialslopes.out")

  open(unit=24, file="plots/MLMCfuncts/spatialslopes.out")

  1130 format(f15.7)
  1131 format("slope for ")
  do ifunct=1,num_ufunct
    write(24,1131,advance="no")
    if(def_ufunct(ifunct,3)==1) then
      write(24,1121) def_ufunct(ifunct,1),def_ufunct(ifunct,2)
    elseif(def_ufunct(ifunct,3)==2) then
      write(24,1122) def_ufunct(ifunct,1),def_ufunct(ifunct,2)
    elseif(def_ufunct(ifunct,3)==3) then
      write(24,1123) def_ufunct(ifunct,1)
    endif
    do igap=1,spatial_Level-1
      do ibase=0,spatial_Level-igap-1
        slope = spatiallogslope(err_ufunctional,ifunct,ibase,ibase+igap)
        write(24,1130,advance="no") slope
      enddo
      write(24,1124)
    enddo
  enddo

  close(unit=24)      


  !Here we print the actual functional error data so that it can be examined
  !by hand and/or plotted.
  call system("test -e plots/MLMCfuncts/spatialconv.out && rm plots/MLMCfuncts/spatialconv.out")

  open(unit=24, file="plots/MLMCfuncts/spatialconv.out")
  1120 format("#ilevel     num of cells   ")
  1121 format("L1 cell",i5," to",i5,"  ")
  1122 format("L2 cell",i5," to",i5,"  ")
  1123 format(" center cell",i5,"     ")
  1124 format(" ")
  write(24,1120,advance="no")
  do ifunct=1,num_ufunct
    if(def_ufunct(ifunct,3)==1) then
      write(24,1121,advance="no") def_ufunct(ifunct,1),def_ufunct(ifunct,2)
    elseif(def_ufunct(ifunct,3)==2) then
      write(24,1122,advance="no") def_ufunct(ifunct,1),def_ufunct(ifunct,2)
    elseif(def_ufunct(ifunct,3)==3) then
      write(24,1123,advance="no") def_ufunct(ifunct,1)
    endif
  enddo
  write(24,1124)

  1125 format(i5,"     ",i13,"    ")
  1126 format(es14.7,"        ")
  do ilevel=0,spatial_Level-1
    write(24,1125,advance="no") ilevel,numcellsLevel0*nextLevelFactor**ilevel
    do ifunct=1,num_ufunct
      write(24,1126,advance="no") err_ufunctional(ifunct,1,ilevel)
    enddo
    write(24,1124)
  enddo
  close(unit=24)    

  end subroutine spatial_calcerr_print


  function spatiallogslope(err_ufunctional2,ifunct,ilevel1,ilevel2)
  !This function finds the log-log slope of spatial convergence of functionals
  !err_ufunctional has a 2 because it is passed by reference, and the new array
  !is starts at 1, not 0.  The '2' denotes this difference for understandability.
  !Each 'ilevel#' has a plus one to compensate for this offset.
  use MLMCvars, only: numcellsLevel0, nextLevelFactor
  real(8) :: err_ufunctional2(:,:,:)
  real(8) :: spatiallogslope
  integer :: ifunct,ilevel1,ilevel2

  spatiallogslope = ( log(err_ufunctional2(ifunct,1,ilevel1+1))                  - &
                      log(err_ufunctional2(ifunct,1,ilevel2+1))                ) / &
                    ( log(real(numcellsLevel0*nextLevelFactor**(ilevel1+1),8))   - &
                      log(real(numcellsLevel0*nextLevelFactor**(ilevel2+1),8)) )
  end function spatiallogslope



  subroutine UQ_iterconv( icase )
  !This subroutine performs an iterative convergence study for all chosen functionals in an
  !effort to determine the iterative convergence parameter for each problem, solve, and QoI
  use genRealzvars, only: s
  use MLMCvars, only: Q_ufunctional, num_ufunct, ncellwidth, numcellsLevel0
  use FEDiffSn, only: setflvarspassedtrue, a, fliterstudy, max_iter, flnewiter, fliterstudy
  integer :: icase, iiter
  real(8), allocatable :: trarray3(:,:,:)

  if(allocated(Q_ufunctional)) deallocate(Q_ufunctional)
  allocate(Q_ufunctional(num_ufunct,1,0:0))
  Q_ufunctional = 0.0d0
  if(allocated(ncellwidth)) deallocate(ncellwidth)
  allocate(ncellwidth(0:0))
  ncellwidth = 0.0d0
  ncellwidth = s/real(numcellsLevel0,8)
  call setflvarspassedtrue                             !tell FE mod to accept input from here
  a             = s

  max_iter = 0
  flnewiter = .true.
  fliterstudy = .true.
  do while(flnewiter)                       !solve until FEDiffSn solver says it's converged
    call sampleconvInput( 0 )               !samples average values, ilevel=0 so to set # of cells
    call solveSamples( max_iter,1,icase )   !solves QoIs, ilevel=max_iter, isamp=1

    if(flnewiter) then
      !add new max_iter level to Q_ufunctional and initialize it
      call move_alloc(Q_ufunctional,trarray3)
      allocate(Q_ufunctional(size(trarray3(:,1,1)),size(trarray3(1,:,1)),0:size(trarray3(1,1,:))))
      Q_ufunctional = 0.0d0
      Q_ufunctional(:,:,0:size(trarray3(1,1,:))-1) = trarray3
      deallocate(trarray3)
    endif

    max_iter = max_iter + 1
  enddo
  call iter_calcerr_print                   !calc and print errs for functs
  end subroutine UQ_iterconv









  subroutine MLMCallocate
  !Allocate arrays here.  Try to only allocate in order to make paralellization easier later.
  use MLMCvars, only: numMLMCcells, M_optsamps, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, bnumMLMCsamps, numcellsLevel0, ncellwidth, num_ufunct, &
                      MLMCerrest

  if(.not.allocated(numMLMCcells)) allocate(numMLMCcells(0:0))
  if(.not.allocated(M_optsamps)) allocate(M_optsamps(3,0:0))
  if(.not.allocated(MLMCerrest)) allocate(MLMCerrest(num_ufunct))

  if(.not.allocated(Q_ufunctional)) allocate(Q_ufunctional(num_ufunct,bnumMLMCsamps,0:0))
  if(.not.allocated(G_ufunctional)) allocate(G_ufunctional(num_ufunct,bnumMLMCsamps,0:0))
  if(.not.allocated(Gave)) allocate(Gave(num_ufunct,0:0))
  if(.not.allocated(Gvar)) allocate(Gvar(num_ufunct,0:0))
  if(.not.allocated(ncellwidth)) allocate(ncellwidth(0:0))

  end subroutine MLMCallocate


  subroutine MLMCinitialize( flMLMC, Level )
  !Initialize/load values to variables here.
  use genRealzvars, only: s
  use MLMCvars, only: numMLMCcells, numcellsLevel0, MLMC_failprob, C_alpha, &
                      M_optsamps, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, bnumMLMCsamps, ncellwidth, MLMCerrest
  use FEDiffSn, only: a
  use utilities, only: erfi
  integer :: Level
  logical :: flMLMC

  Level  = 0
  flMLMC = .true.  

  C_alpha = sqrt(2.0d0)*erfi(1.0d0-MLMC_failprob)

  numMLMCcells    = numcellsLevel0
  M_optsamps(1,0) = bnumMLMCsamps
  M_optsamps(2,0) = 0
  M_optsamps(3,0) = 0

  MLMCerrest    = 0.0d0

  Q_ufunctional = 0.0d0
  G_ufunctional = 0.0d0
  Gave          = 0.0d0
  Gvar          = 0.0d0

  ncellwidth    = s/real(numcellsLevel0,8)
  a             = s

  end subroutine MLMCinitialize


  subroutine MLMCdeallocate
  !Deallocate arrays here.
  use MLMCvars, only: numMLMCcells, M_optsamps, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, ncellwidth, MLMCerrest

  if(allocated(numMLMCcells)) deallocate(numMLMCcells)
  if(allocated(M_optsamps)) deallocate(M_optsamps)
  if(allocated(MLMCerrest)) deallocate(MLMCerrest)

  if(allocated(Q_ufunctional)) deallocate(Q_ufunctional)
  if(allocated(G_ufunctional)) deallocate(G_ufunctional)
  if(allocated(Gave)) deallocate(Gave)
  if(allocated(Gvar)) deallocate(Gvar)
  if(allocated(ncellwidth)) deallocate(ncellwidth)

  end subroutine MLMCdeallocate





  subroutine MLMCaddLevel( Level )
  !Increase size of arrays to handle new Level and in case of response function
  !the associated number of cells as well.
  use MLMCvars, only: numMLMCcells, nextLevelFactor, Q_ufunctional, &
                      G_ufunctional, Gave, Gvar, M_optsamps, bnumMLMCsamps, &
                      ncellwidth
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
  allocate(M_optsamps(3,0:size(tiarray2(1,:))))
  M_optsamps = 0
  M_optsamps(:,0:size(tiarray2(1,:))-1) = tiarray2
  M_optsamps(1,Level) = bnumMLMCsamps
  deallocate(tiarray2)

  !add new Level to Q_ufunctional and initialize it
  call move_alloc(Q_ufunctional,trarray3)
  allocate(Q_ufunctional(size(trarray3(:,1,1)),size(trarray3(1,:,1)),0:size(trarray3(1,1,:))))
  Q_ufunctional = 0.0d0
  Q_ufunctional(:,:,0:size(trarray3(1,1,:))-1) = trarray3
  deallocate(trarray3)

  !add new Level to G_ufunctional and initialize it
  call move_alloc(G_ufunctional,trarray3)
  allocate(G_ufunctional(size(trarray3(:,1,1)),size(trarray3(1,:,1)),0:size(trarray3(1,1,:))))
  G_ufunctional = 0.0d0
  G_ufunctional(:,:,0:size(trarray3(1,1,:))-1) = trarray3
  deallocate(trarray3)

  !add new Level to Gave and initialize it
  call move_alloc(Gave,trarray2)
  allocate(Gave(size(trarray2(:,1)),0:size(trarray2(1,:))))
  Gave = 0.0d0
  Gave(:,0:size(trarray2(1,:))-1) = trarray2
  deallocate(trarray2)

  !add new Level to Gvar and initialize it
  call move_alloc(Gvar,trarray2)
  allocate(Gvar(size(trarray2(:,1)),0:size(trarray2(1,:))))
  Gvar = 0.0d0
  Gvar(:,0:size(trarray2(1,:))-1) = trarray2
  deallocate(trarray2)

  !add new Level to ncellwidth and initialize it
  call move_alloc(ncellwidth,trarray1)
  allocate(ncellwidth(0:size(trarray1(:))))
  ncellwidth = 0.0d0
  ncellwidth(0:size(trarray1(:))-1) = trarray1
  ncellwidth(Level) = ncellwidth(Level-1) / real(nextLevelFactor,8)
  deallocate(trarray1)

  end subroutine MLMCaddLevel


  subroutine MLMCaddSamples( Level )
  !Increase size of arrays to handle new optimal number of samples.
  use MLMCvars, only: Q_ufunctional, G_ufunctional, M_optsamps
  integer :: Level

  real(8), allocatable :: trarray3(:,:,:)

  !add space for new samples to Q_functional and initialize the space
  call move_alloc(Q_ufunctional,trarray3)
  allocate(Q_ufunctional(size(trarray3(:,1,1)),maxval(M_optsamps),0:size(trarray3(1,1,:))-1))
  Q_ufunctional = 0.0d0
  Q_ufunctional(:,1:size(trarray3(1,:,1)),:) = trarray3
  deallocate(trarray3)

  !add space for new samples to G_ufunctional and initialize the space
  call move_alloc(G_ufunctional,trarray3)
  allocate(G_ufunctional(size(trarray3(:,1,1)),maxval(M_optsamps),0:size(trarray3(1,1,:))-1))
  G_ufunctional = 0.0d0
  G_ufunctional(:,1:size(trarray3(1,:,1)),:) = trarray3
  deallocate(trarray3)

  end subroutine MLMCaddSamples



  !3 Using V~s(estimated variance), compute optimal M~s(estimated opt num of samples)
  subroutine MLMCcomputeOptSamps( Level )
  use MLMCvars, only: MLMC_TOLsplit, MLMC_TOL, C_alpha, Gave, Gvar, M_optsamps, &
                      ncellwidth,linsolveEff,numDimensions, def_ufunct
  integer :: Level

  integer :: ilevel,ifunct
  real(8) :: workterm
  real(8) :: accterm !last term in estimate, sum over levels, sum as go

  workterm = 0.0d0
  accterm  = 0.0d0

  do ifunct=1,size(def_ufunct(:,1)) !decide for which functional to set opt num of samps
    if(def_ufunct(ifunct,4)==1) exit
  enddo

  do ilevel=0,Level
    !solve for work term
    workterm = ncellwidth(ilevel)**(-linsolveEff*numDimensions)
    !accumlate next part of accumulating term
    accterm = accterm + sqrt(Gvar(ifunct,ilevel)*workterm)
    !calculate new optimal samples estimate
    M_optsamps(3,ilevel) = ceiling(  (MLMC_TOLsplit*MLMC_TOL/C_alpha)**(-2.0d0) * &
                                     sqrt(abs(Gvar(ifunct,ilevel)/workterm)) * &
                                     accterm                             )                         

    !use larger of new value and old value
    M_optsamps(1,ilevel) = max(M_optsamps(2,ilevel),M_optsamps(3,ilevel))
    !if new level's estimate higher, make backwards compatible
    if(ilevel/=0) then
      if(M_optsamps(1,ilevel) > M_optsamps(1,ilevel-1)) &
         M_optsamps(1,0:ilevel-1) = M_optsamps(1,ilevel)
    endif

    print *,"theoretical optimal M for level",ilevel,":",M_optsamps(3,ilevel)
    print *,"practical   optimal M for level",ilevel,":",M_optsamps(1,ilevel)
  enddo

  end subroutine MLMCcomputeOptSamps




  !4 Evaluate any new samples needed
  subroutine MLMCevalNewSamps( Level,icase )
  !Evaluate any new samples needed, recompute functional values
  use MLMCvars, only: M_optsamps, Gave, Gvar, G_ufunctional, def_ufunct
  use rngvars, only: setrngappnum, rngappnum, rngstride
  use utilities, only: mean_and_var_p
  use mcnp_random, only: RN_init_particle
  integer :: Level, icase

  integer :: ifunct,ilevel,isamp, isamplow

  do ilevel = 0,Level                                          !search each level
    if( M_optsamps(1,ilevel)>M_optsamps(2,ilevel) ) then       !if new samps to compute
      isamplow = max(M_optsamps(2,ilevel)+1,1)                 !set lowest samp num
print *,"isamplow:",isamplow
      do isamp = isamplow,M_optsamps(1,ilevel)                 !cycle through samps to compute
        call setrngappnum('MLMCsamp')                          !set rng unique to sample
        call RN_init_particle( int(rngappnum*rngstride+isamp,8) )
        call sampleInput( ilevel )                             !update solver input info
        call solveSamples( ilevel,isamp,icase )                !solves QoIs
      enddo
    endif
    do ifunct=1,size(def_ufunct(:,1))
      call mean_and_var_p( G_ufunctional(ifunct,:,ilevel),&             !solve ave and var of functionals
                           size(G_ufunctional(ifunct,:,ilevel)),Gave(ifunct,ilevel),Gvar(ifunct,ilevel) )
    enddo
    M_optsamps(2,ilevel) = M_optsamps(1,ilevel)                !save old # of opt samps
  enddo
  print *,"Number of samples solved :",M_optsamps(1,:)
  print *,"Theoretical ideal # samps:",M_optsamps(3,:)
  end subroutine MLMCevalNewSamps



  subroutine sampleInput( ilevel )
  !This subroutine samples input parameters for this solve and passes them as needed
  use genSampvars, only: specialprob, nummat, param1, param2, param1_mean, param1_uncert, &
                         param2_mean, param2_uncert
  use MLMCvars, only: numMLMCcells
  use FEDiffSn, only: sigt, c, numcells,      flvarspassed
  use mcnp_random, only: rang

  integer :: ilevel
  real(8) :: loc_sigs, loc_siga
  logical :: flsiga            !have I already sampled siga?

  !The following preparation of sigt and c is designed for use with FEDiffSn, 
  !and is written on the assumption that there is only 1 material present.
  !When similar parsing is written later, code with that in mind.
  flsiga = .false.

  !sample and prepare sigt
  if( param1(1)=='sigt' ) then
    if( param1(2)=='sigt1-abs' ) then
      sigt = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_uncert(1)
    elseif( param1(2)=='sigt1-frac' ) then
      sigt = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_mean(1)*param1_uncert(1)
    endif
  elseif( param1(1)=='sigs' ) then
    if( param1(2)=='sigs1-abs' ) then
      loc_sigs = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_uncert(1)
    elseif( param1(2)=='sigs1-frac' ) then
      loc_sigs = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_mean(1)*param1_uncert(1)
    endif
    if( param2(2)=='siga1-abs' ) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='siga1-frac' ) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
    flsiga = .true.
    sigt = loc_sigs + loc_siga
  endif

  !sample and prepare c
  if( param2(1)=='c' ) then
    if( param2(2)=='c1-abs' ) then
      c = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='c1-frac' ) then
      c = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
  elseif( param2(1)=='siga' ) then
    if( param2(2)=='siga1-abs' .and. .not.flsiga) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='siga1-frac' .and. .not.flsiga) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
    c = merge(1.0d0-loc_siga/sigt,loc_sigs/(loc_sigs+loc_siga),.not.flsiga)
  elseif( param2(1)=='sigs' ) then
    if( param2(2)=='sigs1-abs' ) then
      loc_sigs = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='sigs1-frac' ) then
      loc_sigs = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
    c = loc_sigs/sigt
  endif

  if( c>1.0d0 ) then
    stop "--'c' sampled as greater than 1, check your input!"
  endif

  numcells = numMLMCcells(ilevel)

  end subroutine sampleInput



  subroutine sampleconvInput( ilevel )
  !Used in spatial convergence studies and iterative convergence studies, 
  !this subroutine samples average input parameters and passes them off
  use genSampvars, only: specialprob, nummat, param1, param2, param1_mean, &
                         param2_mean
  use MLMCvars, only: nextLevelFactor, numcellsLevel0
  use FEDiffSn, only: sigt, c, numcells,      flvarspassed

  integer :: ilevel
  real(8) :: loc_sigs, loc_siga
  logical :: flsiga            !have I already sampled siga?

  !The following preparation of sigt and c is designed for use with FEDiffSn, 
  !and is written on the assumption that there is only 1 material present.
  !When similar parsing is written later, code with that in mind.
  flsiga = .false.

  !sample and prepare sigt
  if( param1(1)=='sigt' ) then
    if( param1(2)=='sigt1-abs' ) then
      sigt = param1_mean(1)
    elseif( param1(2)=='sigt1-frac' ) then
      sigt = param1_mean(1)
    endif
  elseif( param1(1)=='sigs' ) then
    if( param1(2)=='sigs1-abs' ) then
      loc_sigs = param1_mean(1)
    elseif( param1(2)=='sigs1-frac' ) then
      loc_sigs = param1_mean(1)
    endif
    if( param2(2)=='siga1-abs' ) then
      loc_siga = param2_mean(1)
    elseif( param2(2)=='siga1-frac' ) then
      loc_siga = param2_mean(1)
    endif
    flsiga = .true.
    sigt = loc_sigs + loc_siga
  endif

  !sample and prepare c
  if( param2(1)=='c' ) then
    if( param2(2)=='c1-abs' ) then
      c = param2_mean(1)
    elseif( param2(2)=='c1-frac' ) then
      c = param2_mean(1)
    endif
  elseif( param2(1)=='siga' ) then
    if( param2(2)=='siga1-abs' .and. .not.flsiga) then
      loc_siga = param2_mean(1)
    elseif( param2(2)=='siga1-frac' .and. .not.flsiga) then
      loc_siga = param2_mean(1)
    endif
    c = merge(1.0d0-loc_siga/sigt,loc_sigs/(loc_sigs+loc_siga),.not.flsiga)
  elseif( param2(1)=='sigs' ) then
    if( param2(2)=='sigs1-abs' ) then
      loc_sigs = param2_mean(1)
    elseif( param2(2)=='sigs1-frac' ) then
      loc_sigs = param2_mean(1)
    endif
    c = loc_sigs/sigt
  endif

  if( c>1.0d0 ) then
    stop "--'c' sampled as greater than 1, check your input!"
  endif

  numcells = numcellsLevel0*nextLevelFactor**ilevel

  end subroutine sampleconvInput




  subroutine solveSamples( ilevel,isamp,icase )
  !This subroutine drives the solver for a new response function flux for a new sample,
  !collects all desired Quantities of Interest from that solution, 
  !and deallocates module variables from FEDiffSn.
  !For the iter conv study, 'ilevel' here is max_iter, and functional ilevel is 0.
  use MLMCvars, only: Q_ufunctional, G_ufunctional, nextLevelFactor, numcellsLevel0, &
                      ncellwidth, def_ufunct, detMLMC, MLMCcases
  use FEDiffSn, only: FEmain,FEDiffSn_externaldeallocate, &
                      solve,phidiff,phiSnl,phiSnr,phiDSAl,phiDSAr

  integer :: icase, ifunct, isamp, ilevel,  icell, firstcell, lastcell, middlecell, iilevel
  real(8) :: cellwidth
  real(8), allocatable :: flux(:)

  call FEmain

  !collect flux from FEmain results
  if(allocated(flux)) deallocate(flux)
  if(solve(3)==1) then
    allocate(flux(size(phiDSAl)))
    flux = 0.0d0
    flux = (phiDSAl + phiDSAr) / 2.0d0
  elseif(solve(2)==1) then
    allocate(flux(size(phiSnl)))
    flux = 0.0d0
    flux = (phiSnl + phiSnr) / 2.0d0
  elseif(solve(1)==1) then
    allocate(flux(size(phidiff)))
    flux = 0.0d0
    flux = ( phidiff(1:size(phidiff)-1)+phidiff(2:size(phidiff)) )/2.0d0
  endif

  do ifunct=1,size(def_ufunct(:,1))  !solve Q_ufunctional for each specified functional

    iilevel    =  merge(0,ilevel,MLMCcases(icase)=='iter') !to 0 for iter conv study (ilevel=max_iter)
    firstcell  = (def_ufunct(ifunct,1)-1)*nextLevelFactor**iilevel+1
    lastcell   =  def_ufunct(ifunct,2)   *nextLevelFactor**iilevel
    middlecell =  firstcell - 1 +(nextLevelFactor**iilevel+1)/2
    cellwidth  =  ncellwidth(iilevel)


    select case (def_ufunct(ifunct,3))

      case (1) !L1 norm
        Q_ufunctional(ifunct,isamp,ilevel) = 0.0d0
        do icell=firstcell,lastcell
          Q_ufunctional(ifunct,isamp,ilevel) = Q_ufunctional(ifunct,isamp,ilevel) + &
                                                flux(icell)   *cellwidth
        enddo

      case (2) !L2 norm
        Q_ufunctional(ifunct,isamp,ilevel) = 0.0d0
        do icell=firstcell,lastcell
          Q_ufunctional(ifunct,isamp,ilevel) = Q_ufunctional(ifunct,isamp,ilevel) + &
                                                flux(icell)**2*cellwidth
        enddo
        Q_ufunctional(ifunct,isamp,ilevel) = sqrt(Q_ufunctional(ifunct,isamp,ilevel))

      case (3) !center value
        Q_ufunctional(ifunct,isamp,ilevel) = flux(firstcell)
    end select

    !solve G_ufunctional if performing deterministic MLMC
    if( detMLMC=='detMLMC' ) then
      if( ilevel==0 ) then
        G_ufunctional(ifunct,isamp,ilevel) = Q_ufunctional(ifunct,isamp,ilevel)
      else
        G_ufunctional(ifunct,isamp,ilevel) = Q_ufunctional(ifunct,isamp,ilevel) - &
                                             Q_ufunctional(ifunct,isamp,ilevel-1)
      endif

      if(mod(isamp,1000)==0) print *,"G_ufunctional(",ifunct,",",isamp,",",ilevel,"):",&
                                      G_ufunctional(ifunct,isamp,ilevel)
    endif

  enddo

  call FEDiffSn_externaldeallocate

  end subroutine solveSamples



  !5 test total error, set flMLMC==.false.?
  function MLMCcalcErrEst( Level,ifunct )
  use MLMCvars, only: Gave, Gvar, C_alpha, ncellwidth, nextLevelFactor, M_optsamps
             
  real(8) :: MLMCcalcErrEst, err1, C_w, err2, Vsum
  integer :: Level, ilevel, ifunct
  real(8) :: qq = 2.0d0  !Lect 10, pg 3, >=1, order of error, example of =2

  C_w  = max( abs(Gave(ifunct,Level  ))/(ncellwidth(Level  )**qq*(nextLevelFactor**qq-1)), &
              abs(Gave(ifunct,Level-1))/(ncellwidth(Level-1)**qq*(nextLevelFactor**qq-1))    )
  err1 = C_w * ncellwidth(Level)**qq

  Vsum = 0.0d0
  do ilevel=0,Level
    Vsum = Vsum + abs(Gvar(ifunct,ilevel)/M_optsamps(1,ilevel))
  enddo
  err2 = C_alpha * sqrt(Vsum)

  MLMCcalcErrEst = err1 + err2

  print *,"functional:",ifunct
  print *,"  err1, disc err:",err1
  print *,"  err2, MC   err:",err2
  end function MLMCcalcErrEst



  subroutine MLMCprintfunctionaldata
  !Prints data in file for each type of functional chosen
  use genRealzvars, only: s
  use MLMCvars, only: Gave, Gvar, numcellsLevel0, def_ufunct, MLMCerrest, MLMC_TOL

  integer :: L1s, L2s, centers, ifunct

  !count number of each functional to be printed
  L1s     = 0
  L2s     = 0
  centers = 0
  do ifunct=1,size(def_ufunct(:,1))
    if(def_ufunct(ifunct,3)==1) L1s     = L1s + 1
    if(def_ufunct(ifunct,3)==2) L2s     = L2s + 1
    if(def_ufunct(ifunct,3)==3) centers = centers + 1
  enddo

  1110 format("    converged ")
  1111 format("    not converged ")

  !print file for L1s
  call system("test -e plots/MLMCfuncts/MLMCL1s.out && rm plots/MLMCfuncts/MLMCL1s.out")
  if(L1s>0) then
    1100 format("#starting cell, ending cell,  cell center, L1 of flux, estimated error, err tol")
    1101 format(i5,i5,f15.7,f15.7,f15.7,f15.7)
    open(unit=24, file="plots/MLMCfuncts/MLMCL1s.out")
    write(24,1100)
    do ifunct=1,size(def_ufunct(:,1))
      if(def_ufunct(ifunct,3)==1) then
        write(24,1101,advance="no") def_ufunct(ifunct,1),def_ufunct(ifunct,2),&
             ( real(def_ufunct(ifunct,1)+def_ufunct(ifunct,2),8) -1.0d0)*s*0.5d0/real(numcellsLevel0,8), &
             sum(Gave(ifunct,:)),MLMCerrest(ifunct),MLMC_TOL
        if(MLMCerrest(ifunct)<MLMC_TOL) then
          write(24,1110)
        else
          write(24,1111)
        endif
      endif
    enddo
    close(unit=24)    
  endif

  !print file for L2s
  call system("test -e plots/MLMCfuncts/MLMCL2s.out && rm plots/MLMCfuncts/MLMCL2s.out")
  if(L2s>0) then
    1102 format("#starting cell, ending cell,  cell center, L2 of flux, estimated error, err tol")
    1103 format(i5,i5,f15.7,f15.7,f15.7,f15.7)
    open(unit=24, file="plots/MLMCfuncts/MLMCL2s.out")
    write(24,1102)
    do ifunct=1,size(def_ufunct(:,1))
      if(def_ufunct(ifunct,3)==2) then
        write(24,1103,advance="no") def_ufunct(ifunct,1),def_ufunct(ifunct,2),&
             ( real(def_ufunct(ifunct,1)+def_ufunct(ifunct,2),8) -1.0d0)*s*0.5d0/real(numcellsLevel0,8), &
             sum(Gave(ifunct,:)),MLMCerrest(ifunct),MLMC_TOL
        if(MLMCerrest(ifunct)<MLMC_TOL) then
          write(24,1110)
        else
          write(24,1111)
        endif
      endif
    enddo
    close(unit=24)    
  endif

  !print file for centers
  call system("test -e plots/MLMCfuncts/MLMCcenters.out && rm plots/MLMCfuncts/MLMCcenters.out")
  if(centers>0) then
    1104 format("#cell center,    flux at point,     estimated error, err tol")
    1105 format(f15.7,f15.7,f15.7,f15.7)
    open(unit=24, file="plots/MLMCfuncts/MLMCcenters.out")
    write(24,1104)
    do ifunct=1,size(def_ufunct(:,1))
      if(def_ufunct(ifunct,3)==3) then
        write(24,1105,advance="no") ( real(def_ufunct(ifunct,1),8) -0.5d0)*s/real(numcellsLevel0,8), &
                       sum(Gave(ifunct,:)),MLMCerrest(ifunct),MLMC_TOL
        if(MLMCerrest(ifunct)<MLMC_TOL) then
          write(24,1110)
        else
          write(24,1111)
        endif
      endif
    enddo
    close(unit=24)    
  endif

  end subroutine

end module MLMC

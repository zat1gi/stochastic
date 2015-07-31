program stochastic
  use mcnp_random
  use timeman
  use utilities
  use Loadcase
  use genRealz
  use radtransMC
  use KLresearch
  use KLreconstruct
  use KLmeanadjust
  use MLMC
  use FEDiffSn

  use genRealzvars
  use timevars, only: t1
  use KLvars, only: KLrnumRealz, KLrprintat, KLres, KLrec, KLnoise, KLadjust
  use MCvars, only: pltflux, radMC, radWood, KLWood, MCcaseson, probtype
  use MLMCvars, only: MLMCcaseson, MLMCcases

  implicit none
  ! pass by reference
  integer :: j,icase !current realization, current MCtransport case

  real(8) :: rand1,rand2,randGauss,randGauss2(2)

  !test routine for new subprograms OneGaussrandnum and TwoGaussrandnums
  do j=1,10
  rand1 = rang()
  rand2 = rang()
  randGauss = OneGaussrandnum(rand1,rand2)
  call TwoGaussrandnums(rand1,rand2,randGauss2)
  print *,"randGauss : ",randGauss
  print *,"randGauss2: ",randGauss2
  enddo
  stop

  !!read parameters
  call cpu_time(t1)
  call readinputstoc
  call Acase_load !need to load these to test
  call testinputstoc

  !!allocate/prepare global parameters
  call global_allocate
  call Acase_print


  !!Perform KL research
  if(KLres=='yes') then
    call KL_eigenvalue      !solves for eigenmode vars
    call KL_Correlation     !calcs and can plot spacial correlation funcs
    call KL_collect         !collects xi values over realizations
    call genReal_stats      !performs stats on above realizations
    call KL_Cochart         !creates plots of variance kept to total variance
    call KL_eval            !creates xi distributions from xi values
    if(KLnoise=='yes') call KL_Noise !does something with xi distributions
    call reset_genRealtals  !resets Markov realz stats for next round of creation
  endif

  !!Perform KL reconstructions
  if(KLrec=='yes') then
    call KLrmeshgen         !creates mesh for fixed x and xi material constructions
    call KLrgenrealz        !selects array of random variables xi
    if(KLadjust=='yes') call KLadjustmean !adjusts mean after lopping neg cross sections
    call KLrplotrealz       !plots reconstructed realiztions
  endif



  !!Perform UQ-MC for transport problems  
  if( sum(MCcaseson)>0 .and. probtype=='material') then !perform if at least one case chosen
    do icase=1,size(MCcaseson)       !cycle through possible cases
      if( MCcaseson(icase)==1 ) then !run case if chosen
        call UQ_MC( icase )          !perform transport
      endif
    enddo
  endif
  call MCfluxPrint
  call MCfluxPlot
  call MCLeakage_pdfplot

  !!Perform MLMC with deterministic transport
  if( sum(MLMCcaseson)>0 .and. probtype=='coeffs') then                  !perform if at least one case chosen
    do icase = 1,size(MLMCcaseson)                                       !cycle through possible cases
      if( MLMCcaseson(icase)==1 .and. MLMCcases(icase)=='detMLMC' ) then !run 'detMLMC' if chosen
        call UQ_MLMC( icase )                                            !perform UQ problem
      endif
      if( MLMCcaseson(icase)==1 .and. MLMCcases(icase)=='spatial' ) then !run 'spatial'
        call UQ_spatialconv( icase )                                     !perform spatial convergence study
      endif
      if( MLMCcaseson(icase)==1 .and. MLMCcases(icase)=='iter' ) then    !run 'iter'
        call UQ_iterconv( icase )                                        !perform iter convergence study
      endif
      if( MLMCcaseson(icase)==1 .and. MLMCcases(icase)=='bench' ) then   !run 'bench'
        call UQ_benchmark( icase )                                       !generate benchmark values (conv study)
      endif
    enddo
  endif

  !!print final reports
  if(probtype=='material') then
    call clearreports
    call Acase_print
    call Woodnegstats
    if(sum(MCcaseson)/=0) call MCprintstats
    call timereport
    call finalreport
  endif

end program stochastic

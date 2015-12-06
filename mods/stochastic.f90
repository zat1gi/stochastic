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

  use genRealzvars
  use timevars, only: t1
  use KLvars, only: KLres, KLrec, KLnoise, &
                    Corropts, pltCo
  use MCvars, only: pltflux

  implicit none
  ! pass by reference
  integer :: j !current realization

  !!read parameters
  call cpu_time(t1)
  call read_test_inputstoc

  !!allocate/prepare global parameters
  call global_allocate
  call Acase_print
  call clearreports

  !!Perform KL research
  if(KLres=='yes') then
    call KL_eigenvalue      !solves for eigenmode vars
    if(Corropts(1).ne."noplot" .and. .not.flGBgeom) call KL_Correlation !calc & plot spacial correlation funcs
    call KL_collect         !collects xi values over realizations
    call genReal_stats      !performs stats on above realizations
    if(pltCo(1).ne.'noplot' .and. .not.flGBgeom)    call KL_Cochart !creates plots of var kept to tot var
    call KL_eval            !creates xi distributions from xi values
    if(KLnoise=='yes') call KL_Noise !does something with xi distributions
    call reset_genRealtals  !resets Markov realz stats for next round of creation
  endif

  !!Perform KL reconstructions if no transport to use them
  if(KLrec=='yes' .and. .not.chTrantype=='KLWood' .and. .not.chTrantype=='GaussKL') then
    call KLreconstructions
  endif

  !!Perform UQ-MC for transport problems  
  call UQ_MC          !perform transport
  call MCfluxPrint
  call MCfluxPlot
  call MCLeakage_pdfplot

  !!print final reports
  call Acase_print
  call MCprintstats
  call timereport
  call finalreport

end program stochastic

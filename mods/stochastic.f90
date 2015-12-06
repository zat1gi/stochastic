program stochastic
  use mcnp_random
  use utilities
  use Loadcase
  use genRealz
  use radtransMC
  use KLresearch
  use KLconstruct
  use KLmeanadjust

  use genRealzvars
  use KLvars, only: KLres, KLrec, Corropts, pltCo
  use MCvars, only: pltflux

  implicit none
  ! pass by reference
  integer :: j !current realization

  !!read parameters
  call read_test_inputstoc

  !!allocate/prepare global parameters
  call global_allocate
  call Acase_print
  call clearreports

  !!Perform KL research
  if(KLres=='yes') then
    call KL_eigenvalue      !solves for eigenmode vars
    if(Corropts(1).ne."noplot") call KL_Correlation !calc & plot spacial correlation funcs
    call KL_collect         !collects xi values over realizations
    call genReal_stats      !performs stats on above realizations
    if(pltCo(1).ne.'noplot')    call KL_Cochart !creates plots of var kept to tot var
    call KL_binrandvarvals  !creates xi distributions from xi values
    call reset_genRealtals  !resets Markov realz stats for next round of creation
  endif

  !!Perform KL reconstructions if no transport to use them
  if(KLrec=='yes' .and. .not.chTrantype=='KLWood' .and. .not.chTrantype=='GaussKL') then
    call KLconstructions
  endif

  !!Perform UQ-MC for transport problems  
  call UQ_MC          !perform transport
  call MCfluxPrint
  call MCfluxPlot
  call MCLeakage_pdfplot

  !!print final reports
  call Acase_print
  call MCprintstats
  call finalreport

end program stochastic

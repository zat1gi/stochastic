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
  use KLvars, only: Corropts, pltCo, flLNxscheck, flmeanadjust, pltKLrealz
  use MCvars, only: pltflux, binplot, flfluxplot, chTrantype
  implicit none

  !!read parameters
  call read_test_inputstoc

  !!allocate/prepare global parameters
  call global_allocate
  if(chgeomtype=='binary' .and. Adamscase/=0) call Acase_print
  call clearreports

  !!Perform basic KL options
  if(chTrantype=='KLWood' .or. chTrantype=='GaussKL' .or. &
     Corropts(1).ne.'noplot' .or. pltCo(1).ne.'noplot' .or. pltKLrealz(1).ne.'noplot') then
    call KL_eigenvalue                                 !solves for eigenmode vars
    if(Corropts(1).ne.'noplot'    ) call KL_Correlation!calc & plot spacial correlation funcs
    if(chgeomtype=='binary') call KL_collect           !collects xi values over realizations
    if(pltCo(1).ne.'noplot'       ) call KL_Cochart    !creates plots of var kept to tot var
    if(chgeomtype=='binary') then
      call genReal_stats     !performs stats on above realizations
      call KL_binrandvarvals !creates xi distributions from xi values
      call reset_genRealtals !resets Markov realz stats for next round of creation
    endif
  endif

  !!Enable construction of KL realizations
  if(chTrantype=='KLWood' .or. chTrantype=='GaussKL' .or. pltKLrealz(1).ne.'noplot') then
    call KLrmeshgen   !creates mesh for fixed x and xi material constructions
    call KLrgenrealz  !selects array of random variables xi and tests for negativity
    if(flmeanadjust) call KLadjustmean('scatter') !adjusts scat mean after lopping neg xss
    if(flmeanadjust) call KLadjustmean('absorb')  !adjusts abs mean after lopping neg xss
    if(pltKLrealz(1).ne.'noplot') call KLrplotrealz  !plots reconstructed realizations
    if(chTrantype=='GaussKL' .and. flLNxscheck) call LNxsvalstest  !tests local reconstructin variance
  endif

  !!Perform UQ-MC for transport problems  
  if(.not.chTrantype=='None') then
    call UQ_MC          !perform transport
    if(flfluxplot) call MCfluxPrint
    if(flfluxplot) call MCfluxPlot
    if(.not.binplot=='noplot') call MCLeakage_pdfplot
    !!print final reports
    call MCprintstats
    if(chTrantype=='KLWood' .or. chTrantype=='GaussKL') call finalreport
  endif

end program stochastic

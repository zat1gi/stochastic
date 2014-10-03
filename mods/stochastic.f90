program stochastic
  use genRealzvars
  use mcnp_random
  use timeman
  use utilities
  use Loadcase
  use genRealz
  use radtransMC
  use KLresearch
  use KLreconstruct
  use KLmeanadjust
  use Woodcock

  use timevars, only: t1
  use KLvars, only: KLrnumRealz, KLrprintat, KLres, KLrec, KLnoise
  use MCvars, only: pltflux, Wood, radMC, radWood, KLWood, MCcaseson, &
                    numPosMCmeths, &
                    fluxfaces, fluxall, mu, stocMC_fluxall !for test
  use genRealzvars, only: matLength !for test
  implicit none
  ! pass by reference
  integer :: j,icase !current realization, current MCtransport case
  integer :: seed    !random number seed for overall problem, used once.

  !!read parameters
  call cpu_time(t1)
  call readinputstoc( seed )
  call Acase_load !need to load these to test
  call testinputstoc

  !!allocate/prepare global parameters
  call global_allocate( seed )
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
  if( sum(MCcaseson)>0 ) then        !perform if at least one case chosen
    do icase=1,numPosMCmeths         !cycle through possible cases
      if( MCcaseson(icase)==1 ) then !run case if chosen
        call UQ_MC( icase )          !perform transport
      endif
    enddo
  endif
  call MCfluxPrint
  call MCfluxPlot

!  call MCLeakage_pdfplot !bin for radMC,radWood, needs to be updated


  !!print final reports
  call clearreports
  call Acase_print
  call Woodnegstats
  if(sum(MCcaseson)/=0) call MCprintstats
  call timereport
  call finalreport

stop





  !!genRealz, KLresearch, radtrans, radWood
  if(radWood=='yes') Wood='rad'
  do j=1,numRealz
    call genReal( j,'binary ' )
    if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_collect( j )
    if(radMC=='yes') call radtrans_MCsim( j )
    if(radWood=='yes') call WoodcockMC( j )
  enddo
  call genReal_stats
  if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_stats_plot

  !!radKL transport
  if(KLWood=='yes') then
    do j=1,KLrnumRealz !for Woodcockreconstruct later
      call WoodcockMC( j )
    enddo
  endif


end program stochastic

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
  use MCvars, only: pltflux, allowneg, Wood, radMC, radWood, KLWood, MCcaseson, &
                    numPosMCmeths
  implicit none
  ! pass by reference
  integer :: j,icase !current realization, current MCtransport case
  integer :: seed    !random number seed for overall problem, used once.

  !!read parameters
  call cpu_time(t1)
  call readinputstoc( seed )
  call testinputstoc

  !!allocate/prepare global parameters
  call Acase_load
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
  if( sum(MCcaseson)>0 ) then        !perform if at least one cases chosen
    do icase=1,numPosMCmeths         !cycle through possible cases
      if( MCcaseson(icase)==1 ) then !run case if chosen
        call UQ_MC( icase )          !perform transport
      endif
    enddo
  endif


  !!print final reports
  call Acase_print
  if(allowneg=='yes') call Woodnegstats
  if(radMC=='yes' .or. radWood=='yes' .or. KLWood=='yes') call MCprintstats
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
    if(radMC=='yes' .or. KLres=='yes' .or. radWood=='yes') call radtrans_time( j )
  enddo
  call genReal_stats
  if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_stats_plot

  !!radKL transport
  if(KLWood=='yes') then
    Wood='KL'
    do j=1,KLrnumRealz !for Woodcockreconstruct later
      call WoodcockMC( j )
      if(mod(j,KLrprintat)==0) call KLWood_time( j )
    enddo
  endif
!  if(radMC=='yes' .OR. radWood=='yes' .OR. KLWood=='yes')  call transplot( Adamscase )

  !!MC stats and plots.  Some of this needs to be in UQ wrapper, some just passed to 'MCprintstats'
  if(radMC=='yes') call radtrans_MCoutstats
  if(radWood=='yes') call WoodcockMCoutstats
  if(KLWood=='yes') call WoodcockKLoutstats
  if(pltflux(1)/='noplot') call plot_flux
  call MCLeakage_pdfplot !bin for radMC,radWood


end program stochastic

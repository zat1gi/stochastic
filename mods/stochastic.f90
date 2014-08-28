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

  use timevars, only: t1, runtime
  use KLvars, only: KLrnumRealz, KLrprintat, KLres, KLrec, KLnoise
  use MCvars, only: pltflux, allowneg, Wood, radMC, radWood, KLWood
  implicit none
  ! pass by reference
  integer :: j      !current realization, better to make global?
  integer :: seed   !random number seed for overall problem, used once.
  ! local variables
  real(8) :: seeddum

  !!read parameters
  call cpu_time(t1)
  call readinputstoc( seed )
  call testinputstoc

  !!allocate/prepare global parameters
  call Acase_load
  call Acase_print
  do j=1,seed !advance starting seed
    seeddum = rang()
  enddo

  !!Perform KL research
  if(KLres=='yes')   call KL_eigenvalue
  if(KLres=='yes')   call KL_Correlation
  if(radWood=='yes' .or. KLWood=='yes' .or. radMC=='yes' .or. plotmatdxs/='noplot')&
                           call initialize_fluxplot


  !!Perform Transport with Various UQ Methods
  !loop over different methods that will be used
    !allocate UQ specific variables



  !!genRealz, KLresearch, radtrans, radWood
  if(radWood=='yes') Wood='rad'
  do j=1,numRealz
    call genReal( j )
    if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_collect( j )
    if(radMC=='yes') call radtrans_MCsim( j )
    if(radWood=='yes') call WoodcockMC( j )
    if(KLres=='yes') call KL_collect( j )
    if(radMC=='yes' .or. KLres=='yes' .or. radWood=='yes') call radtrans_time( j )
  enddo
  call genReal_stats
  if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_stats_plot
  if(KLres=='yes') call KL_Cochart
  if(KLres=='yes') call KL_eval
  if(KLnoise=='yes') call KL_Noise

  !!KLreconstructions
  if(KLrec=='yes') then
    call KLrcondition
    do j=1,KLrnumRealz
      call KLrgenrealz( j )
      if(mod(j,KLrprintat)==0) call KLr_time( j )
    enddo
  endif
  if(KLadjust=='yes') call KLadjustmean
  if(KLrec=='yes') call KLreval

  !!radKL transport
  if(KLWood=='yes') then
    Wood='KL'
    do j=1,KLrnumRealz !for Woodcockreconstruct later
      call WoodcockMC( j )
      if(mod(j,KLrprintat)==0) call KLWood_time( j )
    enddo
  endif
!  if(radMC=='yes' .OR. radWood=='yes' .OR. KLWood=='yes')  call transplot( Adamscase )





  !!concluding stats-need to be in UQ wrapper
  call Acase_print
  if(radMC=='yes') call radtrans_MCoutstats
  if(radWood=='yes') call WoodcockMCoutstats
  if(KLWood=='yes') call WoodcockKLoutstats



  !!print final reports
  if(KLWood=='yes' .and. allowneg=='yes') call Woodnegstats
  if(pltflux(1)/='noplot') call plot_flux
  call radtrans_resultplot !bin for radMC,radWood
  if(radMC=='yes' .or. radWood=='yes' .or. KLWood=='yes') call MCprintstats
  call timereport


end program stochastic

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

!These here only for a test
  real(8) :: minpos, maxpos
  integer :: i,loop
  character(12) :: flfluxtallytype
  integer :: ibin
!end test vars

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

  !!print final reports
  call clearreports
  call Acase_print
  call Woodnegstats
  if(sum(MCcaseson)/=0) call MCprintstats
  call timereport
  call finalreport




  !!test loop for MCfluxtally, comparing to older version
  if(allocated(matLength)) deallocate(matLength)
  allocate(matLength(4))
  matLength = (/ 0.0d0, 0.2d0, 0.7d0, 1.0d0 /)
  if(allocated(fluxall)) deallocate(fluxall)
  allocate(fluxall(10,2))
  fluxall = 0.0d0
  if(allocated(fluxfaces)) deallocate(fluxfaces)
  allocate(fluxfaces(11))
  fluxfaces = (/ 0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0 /)
  print *,"fluxall: ",fluxall
  do loop=1,5
    j=1
    flfluxtallytype = 'irrespective'
    if(loop==1) then
      minpos = 0.5d0
      maxpos = 0.71d0
      mu = 0.73d0
    elseif(loop==2) then
      minpos = 0.11d0
      maxpos = 0.99d0
      mu = 0.412d0
    elseif(loop==3) then
      minpos = 0.0d0
      maxpos = 1.0d0
      mu = -0.88d0
    elseif(loop==4) then
      minpos = 0.143d0
      maxpos = 0.205d0
      mu = 0.011d0
    elseif(loop==5) then
      minpos = 0.799d0
      maxpos = 0.95d0
      mu = -0.6d0
    endif
    i=1   !can probably do this part with a ceiling function, once debugged look at again
    do
      if(matLength(i)<=minpos .and. minpos<matLength(i+1)) exit
      i = i+1
    enddo

    call MCfluxtally( j,i,minpos,maxpos,flfluxtallytype )
  enddo
  print *,"fluxall: ",fluxall/(fluxfaces(2)-fluxfaces(1))

  fluxall = fluxall / (fluxfaces(2)-fluxfaces(1))
  fluxall(1,2)  = 1.2d0
  fluxall(2,2)  = 1.2d0
  fluxall(3,2)  = 1.3d0
  fluxall(4,2)  = 1.4d0
  fluxall(5,2)  = 1.2d0
  fluxall(6,2)  = 1.4d0
  fluxall(7,2)  = 1.0d0
  fluxall(8,2)  = 0.9d0
  fluxall(9,2)  = 0.7d0
  fluxall(10,2) = 0.3d0

  if(allocated(stocMC_fluxall)) deallocate(stocMC_fluxall)
  allocate(stocMC_fluxall(10,1,2))
  stocMC_fluxall = 0.0d0

  fluxall = fluxall / 5.0d0
  do ibin=1,10
    call mean_and_var_s( fluxall(ibin,:),2, &
         stocMC_fluxall(ibin,1,1),stocMC_fluxall(ibin,1,2) )
  enddo

  print *,"stocMC_fluxall(:,1,1): ",stocMC_fluxall(:,1,1)
  print *,"stocMC_fluxall(:,1,2): ",stocMC_fluxall(:,1,2)


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

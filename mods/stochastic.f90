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

  use KLvars, only: KLrnumRealz, KLrprintat
  use MCvars, only: pltflux, allowneg, Wood
  implicit none
  real(8) :: runtime,t1,t2,seeddum
  integer :: seed
  character(3) :: KLres,KLrec,radMC,KLnoise,radWood,KLWood
  !--- genRealz variables (new) ---!
  integer :: i,j
  !--- radtransMC variables (new) ---!
  integer      :: o
  real(8),allocatable :: aveRefl(:),devRefl(:),relRefl(:)
  real(8),allocatable :: aveTran(:),devTran(:),relTran(:)
  real(8),allocatable :: aveAbso(:),devAbso(:),relAbso(:)
  real(8),allocatable :: fluxfaces(:),flux(:,:),fflux(:,:),bflux(:,:)
  !--- Woodcock variables (new) ---!
  real(8),allocatable :: Woodf(:,:),radWoodf(:,:),KLWoodf(:,:)
  real(8),allocatable :: fWoodf(:,:),bWoodf(:,:)
  real(8),allocatable :: fradWoodf(:,:),bradWoodf(:,:),fKLWoodf(:,:),bKLWoodf(:,:)

  !!read and prepare parameters
  call cpu_time(t1)
  call readinputstoc(      KLres,KLnoise,&
                           KLrec,radMC,&
                           radWood,KLWood,&
                           seed )

  call testinputstoc(      KLres,KLrec,radWood,&
                           radMC,KLnoise,KLWood )

  call Acase_load
  call Acase_print

  do j=1,seed !advance starting seed
    seeddum = rang()
  enddo

  !!genRealz, KLresearch, radtrans, radWood
  if(KLres=='yes')   call KL_eigenvalue
  if(KLres=='yes')   call KL_Correlation
  if(radWood=='yes' .OR. KLWood=='yes' .OR. radMC=='yes' .or. plotmatdxs/='noplot')&
                           call initialize_fluxplot(&
                           fluxfaces,flux,fflux,bflux,&
                           radMC,radWood,KLWood,radWoodf,KLWoodf,&
                           fradWoodf,bradWoodf,fKLWoodf,bKLWoodf )
  do j=1,numRealz
    call genReal(          j )
    if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_collect( &
                           j,fluxfaces )
    if(radMC=='yes') call radtrans_MCsim( j,&
                           o,fluxfaces,flux,fflux,bflux,s )
    if(radWood=='yes') Wood='rad'
    if(radWood=='yes') call WoodcockMC( j,&
                           fluxfaces,Woodf,radWoodf,KLWoodf,&
                           fWoodf,bWoodf,fradWoodf,bradWoodf,&
                           fKLWoodf,bKLWoodf )
    if(KLres=='yes') call KL_collect( j )
    if(radMC=='yes' .OR. KLres=='yes' .OR. radWood=='yes') call radtrans_time( &
                           radMC,KLres,radWood,j,t1 )
  enddo
  call genReal_stats
  if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_stats_plot( &
                           fluxfaces )
  if(KLres=='yes') call KL_Cochart
  if(KLres=='yes') call KL_eval
  if(KLnoise=='yes') call KL_Noise



  !!KLreconstructions
  if(KLrec=='yes') call KLrcondition
  do j=1,KLrnumRealz
    if(KLrec=='yes') call KLrgenrealz( j,t1 )
    if(mod(j,KLrprintat)==0 .AND. KLrec=='yes') call KLr_time( j,t1)
  enddo
  if(KLadjust=='yes') call KLadjustmean
  if(KLrec=='yes') call KLreval


  !!radKL transport
  do j=1,KLrnumRealz !for Woodcockreconstruct later

    if(KLWood=='yes') Wood='KL'
    if(KLWood=='yes') call WoodcockMC( j,&
                         fluxfaces,Woodf,radWoodf,KLWoodf,&
                         fWoodf,bWoodf,fradWoodf,bradWoodf,&
                         fKLWoodf,bKLWoodf )
    if(mod(j,KLrprintat)==0 .AND. KLWood=='yes') call KLWood_time( j,t1)
  enddo

!  if(radMC=='yes' .OR. radWood=='yes' .OR. KLWood=='yes')  call transplot( Adamscase )



  !!concluding stats
  call Acase_print

  if(radMC=='yes') call radtrans_MCoutstats( flux,fluxfaces,fflux,bflux )
  if(radWood=='yes') call WoodcockMCoutstats( fluxfaces,radWoodf,fradWoodf,bradWoodf )
  if(KLWood=='yes') call WoodcockKLoutstats( fluxfaces,KLWoodf,fKLWoodf,bKLWoodf )
  if(KLWood=='yes' .and. allowneg=='yes') call Woodnegstats
  if(pltflux(1)/='noplot') call plot_flux( radMC,radWood,KLWood )
  call radtrans_resultplot !bin for radMC,radWood

  write(*,*)
  call calc_time_p(        t1,t2,runtime )
  call timereport(         runtime,KLres,KLrec,radMC,radWood,KLWood,KLnoise )

end program stochastic

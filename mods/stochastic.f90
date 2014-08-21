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
  implicit none
  real(8) :: runtime,t1,t2,seeddum
  real(8),allocatable :: time(:) !genR,radMC,radWood,KLnoise,KLcol,KLrec,KLWood
  integer :: seed
  integer, parameter :: ntime = 7
  character(3) :: KLres,KLrec,radMC,KLnoise,radWood,KLWood
  !--- genRealz variables (new) ---!
  integer :: i,j,largesti,numPath(2),pltgenrealznumof,nummatSegs
  real(8) :: P(2),perFirstTally(2),devFirstTally(2),lamc
  real(8) :: matFirstTally(2)=0,sumPath(2),sqrPath(2),avePath(2),devPath(2)
  character(7) :: pltgenrealz(4),plotmatdxs
  integer,allocatable :: matType(:),pltgenrealzwhich(:)
  real(8),allocatable :: matLength(:),matdxs(:,:,:)
  !--- radtransMC variables (new) ---!
  integer      :: numParts,trannprt,o,radtrans_int,pfnumcells
  character(6) :: rodOrplanar,results,sourceType,plotflux(2)
  character(7) :: pltflux(4)
  real(8),allocatable :: reflect(:),transmit(:),absorb(:),initcur(:)
  real(8),allocatable :: aveRefl(:),devRefl(:),relRefl(:)
  real(8),allocatable :: aveTran(:),devTran(:),relTran(:)
  real(8),allocatable :: aveAbso(:),devAbso(:),relAbso(:)
  real(8),allocatable :: fluxfaces(:),flux(:,:),fflux(:,:),bflux(:,:)
  !--- KLresearch variables (new) ---!
  integer      :: numEigs,numSlice,levsrefEig,mostinBin,Corrnumpoints
  real(8)      :: binSmallBound,binLargeBound,sigave,totLength(2),binSize,CoExp
  character(7) :: pltxiBins(4),pltxiBinsgauss,pltEigf(4),pltCo(4),Corropts(2)
  !--- KLreconstruct variables (new) ---!
  integer      :: KLrnumpoints(2),KLrnumRealz,KLrprintat,negcnt
  character(7) :: pltKLrrealz(4)
  integer      :: pltKLrrealznumof
  integer,allocatable :: pltKLrrealzwhich(:,:)
  real(8),allocatable :: KLrx(:),KLrxi(:),pltKLrrealzarray(:,:),KLrxivals(:,:)
  real(8),allocatable :: KLrrandarray(:,:,:),KLrsig(:),KLrxisig(:)
  character(7),allocatable :: pltKLrrealzPointorXi(:)
  !--- Woodcock variables (new) ---!
  integer :: Wood_rej(2),radWood_rej(2),KLWood_rej(2), numpnSamp(2)
  real(8) ::                                 disthold,areapnSamp(4)
  real(8),allocatable :: Woodt(:),   Woodr(:),   Wooda(:)
  real(8),allocatable :: radWoodt(:),radWoodr(:),radWooda(:)
  real(8),allocatable :: KLWoodt(:), KLWoodr(:), KLWooda(:)
  real(8),allocatable :: Woodf(:,:),radWoodf(:,:),KLWoodf(:,:)
  real(8),allocatable :: fWoodf(:,:),bWoodf(:,:)
  real(8),allocatable :: fradWoodf(:,:),bradWoodf(:,:),fKLWoodf(:,:),bKLWoodf(:,:)
  character(3) :: Wood,allowneg,distneg

  !!read and prepare parameters
  call cpu_time(t1)
  call readinputstoc(      numEigs,numSlice,levsrefEig,&
                           binSmallBound,binLargeBound,KLres,KLnoise,&
                           pltxiBins,&
                           pltxiBinsgauss,pltKLrrealzPointorXi,&
                           pltEigf,&
                           KLrnumpoints,KLrnumRealz,KLrprintat,KLrec,&
                           pltKLrrealz,pltKLrrealznumof,pltKLrrealzwhich,&
                           numParts,trannprt,radMC,rodOrplanar,results,&
                           pltgenrealz,pltgenrealznumof,pltgenrealzwhich,&
                           pltCo,plotmatdxs,&
                           Corrnumpoints,Corropts,radWood,KLWood,allowneg,&
                           distneg,plotflux,pfnumcells,pltflux,sourceType,seed )

  call testinputstoc(      numEigs,&
                           KLrnumRealz,KLrprintat,&
                           pltKLrrealznumof,pltKLrrealzwhich,pltEigf,pltxiBins,&
                           pltKLrrealz,trannprt,KLres,KLrec,radWood,&
                           pltgenrealz,pltgenrealznumof,pltgenrealzwhich,&
                           pltCo,radMC,pltKLrrealzPointorXi,&
                           KLrnumpoints,KLnoise,KLWood,pltflux,&
                           sourceType,allowneg,distneg )

  call Acase_load
  call Acase_print(        lamc )

  do j=1,seed !advance starting seed
    seeddum = rang()
  enddo

  !!genRealz, KLresearch, radtrans, radWood
  if(KLres=='yes')   call KL_eigenvalue( numEigs,P,sigave,&
                           levsrefEig,lamc,numSlice,pltEigf,&
                           KLrxivals,KLrnumRealz )
  if(KLres=='yes')   call KL_Correlation( Corropts,Corrnumpoints,numEigs,&
                           lamc,sigave,CoExp,P )
  if(radWood=='yes' .OR. KLWood=='yes' .OR. radMC=='yes' .or. plotmatdxs/='noplot')&
                           call initialize_fluxplot(&
                           pfnumcells,fluxfaces,flux,fflux,bflux,&
                           plotflux,radMC,radWood,KLWood,radWoodf,KLWoodf,&
                           KLrnumRealz,fradWoodf,bradWoodf,fKLWoodf,bKLWoodf )
  do j=1,numRealz
    call genReal(          P,matLength,matType,matFirstTally,&
                           numPath,sumPath,sqrPath,largesti,j,time,ntime,&
                           pltgenrealz,pltgenrealznumof,pltgenrealzwhich,nummatSegs )
    if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_collect( matdxs,&
                           j,matLength,matType,fluxfaces,pfnumcells,nummatSegs )
    if(radMC=='yes') call radtrans_MCsim( j,nummatSegs,numParts,&
                           rodOrplanar,matLength,matType,o,transmit,&
                           reflect,absorb,initcur,time,ntime,radtrans_int,&
                           pfnumcells,fluxfaces,flux,fflux,bflux,plotflux,pltflux,&
                           sourceType,s )
    if(radWood=='yes') Wood='rad'
    if(radWood=='yes') call WoodcockMC( j,matType,matLength,nummatSegs,&
                           time,ntime,numParts,lamc,Wood,numEigs,&
                           radWoodt,radWoodr,radWooda,radWood_rej,KLrnumRealz,&
                           Woodt,Woodr,Wooda,KLWoodt,KLWoodr,KLWooda,Wood_rej,&
                           KLWood_rej,sigave,KLrxivals,rodOrplanar,&
                           fluxfaces,plotflux,pltflux,Woodf,radWoodf,KLWoodf,&
                           pfnumcells,sourceType,fWoodf,bWoodf,fradWoodf,bradWoodf,&
                           fKLWoodf,bKLWoodf,allowneg,numpnSamp,areapnSamp,distneg,&
                           disthold )
    if(KLres=='yes') call KL_collect( nummatSegs,matLength,matType,j,&
                           numEigs,sigave,lamc,totLength,&
                           time,ntime )
    if(radMC=='yes' .OR. KLres=='yes' .OR. radWood=='yes') call radtrans_time( time,&
                           ntime,radMC,KLres,radWood,j,trannprt,t1 )
  enddo
  call genReal_stats(      matFirstTally,perFirstTally,devFirstTally,&
                           P,numPath,sumPath,sqrPath,avePath,devPath,&
                           pltgenrealz,pltgenrealznumof,pltgenrealzwhich )
  if(plotmatdxs/='noplot' .or. pltflux(1)/='noplot') call matdxs_stats_plot( matdxs,&
                           plotmatdxs,fluxfaces,pfnumcells )
  if(KLres=='yes') call KL_Cochart( numEigs,numSlice,P,sigave,lamc,&
                           avePath,totLength,pltCo,&
                           CoExp )
  if(KLres=='yes') call KL_eval( binSmallBound,binLargeBound,&
                           numEigs,&
                           pltxiBins,pltxiBinsgauss,binSize,&
                           mostinBin )
  if(KLnoise=='yes') call KL_Noise( numEigs,&
                           binSmallBound,binLargeBound,binSize,mostinBin,time,ntime )



  !!KLreconstructions
  if(KLrec=='yes') call KLrcondition( KLrx,KLrxi,KLrnumpoints,s )
  do j=1,KLrnumRealz
    if(KLrec=='yes') call KLrgenrealz( sigave,numEigs,lamc,KLrx,&
                           KLrnumpoints,j,KLrnumRealz,&
                           KLrprintat,t1,pltKLrrealz,time,ntime,negcnt,&
                           pltKLrrealznumof,pltKLrrealzwhich,pltKLrrealzarray,&
                           KLrrandarray,KLrsig,KLrxisig,KLrxivals,KLrxi )
    if(mod(j,KLrprintat)==0 .AND. KLrec=='yes') call KLr_time( time,ntime,j,&
                           KLrnumRealz,t1)
  enddo
  if(KLadjust=='yes') call KLadjustmean( sigave,numEigs,lamc,&
                           KLrnumRealz,KLrxivals )
  if(KLrec=='yes') call KLreval( KLrnumpoints,pltKLrrealznumof,pltKLrrealzarray,&
                           pltKLrrealz,KLrrandarray,lamc,&
                           KLrx,numEigs,pltKLrrealzwhich,&
                           KLrsig,sigave,pltKLrrealzPointorXi,KLrxi,KLrxisig,&
                           KLrxivals,negcnt )


  !!radKL transport
  do j=1,KLrnumRealz !for Woodcockreconstruct later

    if(KLWood=='yes') Wood='KL'
    if(KLWood=='yes') call WoodcockMC( j,matType,matLength,nummatSegs,&
                         time,ntime,numParts,lamc,Wood,numEigs,&
                         radWoodt,radWoodr,radWooda,radWood_rej,KLrnumRealz,&
                         Woodt,Woodr,Wooda,KLWoodt,KLWoodr,KLWooda,Wood_rej,&
                         KLWood_rej,sigave,KLrxivals,rodOrplanar,&
                         fluxfaces,plotflux,pltflux,Woodf,radWoodf,KLWoodf,&
                         pfnumcells,sourceType,fWoodf,bWoodf,fradWoodf,bradWoodf,&
                         fKLWoodf,bKLWoodf,allowneg,numpnSamp,areapnSamp,distneg,&
                         disthold )
    if(mod(j,KLrprintat)==0 .AND. KLWood=='yes') call KLWood_time( time,ntime,j,&
                           KLrnumRealz,t1)
  enddo

!  if(radMC=='yes' .OR. radWood=='yes' .OR. KLWood=='yes')  call transplot( &
!                          Adamscase )



  !!concluding stats
  call Acase_print(        lamc )

  if(radMC=='yes') call radtrans_MCoutstats( reflect,transmit,absorb,initcur,&
                           numParts,&
                           results,radtrans_int,rodOrplanar,plotflux,pltflux,&
                           pfnumcells,flux,fluxfaces,fflux,bflux,P )
  if(radWood=='yes') call WoodcockMCoutstats( numParts,radWoodt,radWoodr,&
                           radWooda,radWood_rej,plotflux,pltflux,&
                           pfnumcells,fluxfaces,radWoodf,fradWoodf,bradWoodf,P )
  if(KLWood=='yes') call WoodcockKLoutstats( numParts,KLrnumRealz,KLWoodt,KLWoodr,&
                           KLWooda,KLWood_rej,plotflux,pltflux,&
                           pfnumcells,fluxfaces,KLWoodf,fKLWoodf,bKLWoodf,P,numEigs )
  if(KLWood=='yes' .and. allowneg=='yes') call Woodnegstats( negcnt,&
                           numpnSamp,areapnSamp,distneg )
  if(pltflux(1)/='noplot') call plot_flux( plotflux,pltflux,radMC,radWood,KLWood )
  call radtrans_resultplot( reflect,transmit,radWoodt,radWoodr,KLWoodt,KLWoodr )!bin for radMC,radWood

  write(*,*)
  call calc_time_p(        t1,t2,runtime )
  call timereport(         runtime,time,ntime,KLres,KLrec,radMC,radWood,KLWood,&
                           KLnoise )

end program stochastic

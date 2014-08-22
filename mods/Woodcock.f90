module Woodcock
  use mcnp_random
  use timeman
  use radtransMC
  use utilities
  use KLmeanadjust
  implicit none

CONTAINS
  ! print statements in this module use 600-699


  subroutine WoodcockMC( j,matType,matLength,nummatSegs,&
                         time,ntime,numParts,lamc,Wood,&
                         radWoodt,radWoodr,radWooda,radWood_rej,&
                         Woodt,Woodr,Wooda,KLWoodt,KLWoodr,KLWooda,Wood_rej,&
                         KLWood_rej,rodOrplanar,&
                         fluxfaces,plotflux,pltflux,Woodf,radWoodf,KLWoodf,&
                         pfnumcells,sourceType,fWoodf,bWoodf,fradWoodf,bradWoodf,&
                         fKLWoodf,bKLWoodf,allowneg,numpnSamp,areapnSamp,distneg,&
                         disthold )
  use genRealzvars, only: sig, scatrat, lam, s, numRealz
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave, KLrnumRealz
  integer :: j,nummatSegs,numParts,ntime,matType(:),pfnumcells
  integer :: Wood_rej(2),radWood_rej(2),KLWood_rej(2),numpnSamp(2)
  real(8),allocatable :: Woodt(:),   Woodr(:),   Wooda(:)
  real(8),allocatable :: radWoodt(:),radWoodr(:),radWooda(:)
  real(8),allocatable :: KLWoodt(:), KLWoodr(:), KLWooda(:)
  real(8),allocatable :: Woodf(:,:),fWoodf(:,:),bWoodf(:,:)
  real(8) :: matLength(:),time(:),lamc,areapnSamp(4),disthold
  real(8) :: fluxfaces(:),radWoodf(:,:),KLWoodf(:,:)
  real(8) :: fradWoodf(:,:),bradWoodf(:,:),fKLWoodf(:,:),bKLWoodf(:,:)
  character(3) :: Wood,allowneg,distneg
  character(6) :: rodOrplanar,plotflux(2),sourceType
  character(7) :: pltflux(4)

  integer :: i,o,nbin,k
  real(8) :: binlength,position,mu,ceilsig,woodrat,actsig,tscatrat
  real(8) :: db,dc,tt1,tt2
  real(8),allocatable :: binmaxind(:),binmaxes(:),fbinmax(:),bbinmax(:)

  character(3) :: print='no'

  call cpu_time(tt1)

  if(j==1 .AND. Wood=='rad') then !initialize
    allocate(Woodt(numRealz))
    allocate(Woodr(numRealz))
    allocate(Wooda(numRealz))
    if(plotflux(2)=='tot') allocate(Woodf(numRealz,pfnumcells))
    if(plotflux(2)=='fb') allocate(fWoodf(numRealz,pfnumcells))
    if(plotflux(2)=='fb') allocate(bWoodf(numRealz,pfnumcells))
    Wood_rej =0
    Woodt    =0.0d0
    Woodr    =0.0d0
    Wooda    =0.0d0
    if(plotflux(2)=='tot') Woodf   =0.0d0
    if(plotflux(2)=='fb') fWoodf   =0.0d0
    if(plotflux(2)=='fb') bWoodf   =0.0d0
  endif
  if(j==1 .AND. Wood=='KL') then !initialize
    allocate(Woodt(KLrnumRealz))
    allocate(Woodr(KLrnumRealz))
    allocate(Wooda(KLrnumRealz))
    if(plotflux(2)=='tot') allocate(Woodf(numRealz,pfnumcells))
    if(plotflux(2)=='fb')  allocate(fWoodf(numRealz,pfnumcells))
    if(plotflux(2)=='fb')  allocate(bWoodf(numRealz,pfnumcells))
    Wood_rej =0
    Woodt    =0.0d0
    Woodr    =0.0d0
    Wooda    =0.0d0
    if(plotflux(2)=='tot') Woodf   =0.0d0
    if(plotflux(2)=='fb') fWoodf   =0.0d0
    if(plotflux(2)=='fb') bWoodf   =0.0d0
    numpnSamp =0
    areapnSamp=0.0d0
    disthold  =0.0d0
  endif




  !select local bin maxes
  if(Wood=='rad') nbin = ceiling(s/lamc)
  if(Wood=='KL')  nbin = numEigs

  allocate(binmaxind(nbin+1))
  allocate(binmaxes(nbin))
  allocate(fbinmax(nbin))
  allocate(bbinmax(nbin))
  binmaxind=0.0d0
  binmaxes=0.0d0
  fbinmax=0.0d0
  bbinmax=0.0d0

  binlength=s/nbin
  binmaxind(1)=0.0d0
  do i=2,nbin+1
    binmaxind(i)=binmaxind(i-1)+binlength
    if(print=='yes') print *,"i:",i,"binmaxind(i):",binmaxind(i),"   nbin",nbin
  enddo

  if(Wood=='rad') call radWood_binmaxes(matLength,matType,nummatSegs,binmaxind,binmaxes,nbin,sig)
  if(Wood=='KL')  call KLWood_binmaxes( j,lamc,&
                                        binmaxind,binmaxes,nbin)

  !create forward/backward motion max vectors
  bbinmax(1)=binmaxes(1)
  fbinmax(nbin)=binmaxes(nbin)
  do i=2,nbin
    bbinmax(i) = merge(bbinmax(i-1),binmaxes(i),bbinmax(i-1)>binmaxes(i))
  enddo
  do i=nbin-1,1,-1
    fbinmax(i) = merge(fbinmax(i+1),binmaxes(i),fbinmax(i+1)>binmaxes(i))
  enddo


if(print=='yes') print *,"binmaxind",binmaxind
if(print=='yes') print *,"binmaxes ",binmaxes
if(print=='yes') print *,"fbinmax  ",fbinmax
if(print=='yes') print *,"bbinmax  ",bbinmax
if(print=='yes') print *,
if(print=='yes') print *,
if(print=='yes') print *,

  do o=1,numParts !o-loop, through per particle

    if( sourceType=='left' ) then  !generate source particles
      position=0.0d0
      mu=1.0d0
      if(rodOrplanar=='planar') mu = isoboundmu()
    elseif( sourceType=='intern' ) then
      position= s * rang()
      mu=merge(1.0d0,-1.0d0,rang()>=0.5d0)
      if(rodOrplanar=='planar') mu = newmu()
    endif


    if(print=='yes') print *,
    if(print=='yes') print *,"starting particle ",o
    do    !p-loop, through per particle interaction
      db = merge(s-position,position,mu>=0)/abs(mu) !calc db
      
      ceilsig=merge(ceilsigfunc(binmaxind,fbinmax,position,nbin),& !sel max sig
             ceilsigfunc(binmaxind,bbinmax,position,nbin),mu>=0)

      dc = -log(1-rang())/ceilsig                   !calc dc

      if(print=='yes') print *,"  position   :",position,"    mu     :",mu,&
                               "   db:",db,"   ceilsig:",ceilsig,"   dc:",dc


      if( dc>db ) Wood_rej(1)=Wood_rej(1)+1   !accept path
      if( dc>db .AND. mu>=0 ) then        !tally trans or refl
        Woodt(j) = Woodt(j)+1 !transmit
        if(plotflux(2)=='tot') call adv_pos_col_flux(position,s,fluxfaces,Woodf,&
                                    pfnumcells,plotflux,pltflux,j,mu)
        if(plotflux(2)=='fb') call col_fbflux(position,s,fluxfaces,fWoodf,bWoodf,&
                                   pfnumcells,plotflux,pltflux,j,mu,matType,matLength,nummatSegs)

        if(print=='yes') print *,"                      tally transmit"
        exit
      endif
      if( dc>db .AND. mu<0 ) then
        Woodr(j) = Woodr(j) + 1 !reflect
        if(plotflux(2)=='tot') call adv_pos_col_flux(position,0.0d0,fluxfaces,Woodf,&
                                    pfnumcells,plotflux,pltflux,j,mu)
        if(plotflux(2)=='fb') call col_fbflux(position,0.0d0,fluxfaces,fWoodf,bWoodf,&
                                   pfnumcells,plotflux,pltflux,j,mu,matType,matLength,nummatSegs)
        if(print=='yes') print *,"                      tally reflect"
        exit
      endif

      if(plotflux(2)=='tot') call adv_pos_col_flux(position,position+dc*mu,fluxfaces,Woodf,&
                                  pfnumcells,plotflux,pltflux,j,mu)
      if(plotflux(2)=='fb') call col_fbflux(position,position+dc*mu,fluxfaces,fWoodf,bWoodf,&
                                 pfnumcells,plotflux,pltflux,j,mu,matType,matLength,nummatSegs)
      if(Wood=='rad') woodrat= radWood_actsig(position,matType,matLength,sig)&
                               /ceilsig
      if(Wood=='KL')  woodrat= KLrxi_point(j,lamc,position)&
                               /ceilsig
      if(woodrat>1.0d0) then
        print *,"j: ",j,"  woodrat: ",woodrat
        STOP 'Higher sig sampled in KLWood than ceiling, exiting program'
      endif
      if(Wood=='KL' .and. distneg=='yes' .and. woodrat>0.0d0) then !opt to redist neg xs values
        if(abs(disthold)>woodrat*ceilsig) then
          woodrat = 0.0d0
        else
          woodrat = (woodrat*ceilsig + disthold) / ceilsig
        endif
        disthold = 0.0d0
      endif
      if(woodrat<0.0d0 .and. allowneg=='no') STOP 'Neg number sampled in KLWood, exiting program'
      if(Wood=='KL'  .and. allowneg=='yes') then  !tally data for neg/pos if allowing neg
        if(woodrat<0.0) then
          numpnSamp(2)  =  numpnSamp(2)+1
          areapnSamp(2) = areapnSamp(2)+woodrat*ceilsig          
          if(woodrat*ceilsig<areapnSamp(4)) areapnSamp(4)=woodrat*ceilsig
          if(distneg=='yes') disthold = woodrat*ceilsig
!if(distneg=='yes') print *,"disthold: ",disthold
          woodrat=0.0d0
        else
          numpnSamp(1)  =  numpnSamp(1)+1
          areapnSamp(1) = areapnSamp(1)+woodrat*ceilsig          
          if(woodrat*ceilsig>areapnSamp(3)) areapnSamp(3)=woodrat*ceilsig
        endif
      endif

      if(rang()>woodrat) then                 !not true int?
        Wood_rej(2)=Wood_rej(2)+1 !reject path
        if(print=='yes') print *,"                      reject interaction"
        cycle
      endif

      Wood_rej(1)=Wood_rej(1)+1   !accept path

      if(Wood=='rad') tscatrat=radWood_actscatrat(position,matType,matLength,scatrat)
      if(Wood=='KL')  tscatrat=scatrat(1)
      if(rang()<tscatrat) then !scat or absorb
        if(print=='yes') print *,"                      scatter (cycle)"
        if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
        if(rodOrplanar=='planar') mu = newmu()
        cycle
     else
        Wooda(j) = Wooda(j)+1 !absorb
        if(print=='yes') print *,"                      tally absorb"
        exit
      endif

    enddo  !end p-loop
  enddo  !end o-loop





  if(j==numRealz .AND. Wood=='rad') then  !store away for use in stats
    allocate(radWoodt(numRealz))
    allocate(radWoodr(numRealz))
    allocate(radWooda(numRealz))
    radWood_rej(1)=Wood_rej(1)
    radWood_rej(2)=Wood_rej(2)
    do k=1,numRealz
      radWoodt(k)=Woodt(k)
      radWoodr(k)=Woodr(k)
      radWooda(k)=Wooda(k)
      do i=1,pfnumcells
        if(plotflux(2)=='tot') radWoodf(k,i)=Woodf(k,i)
        if(plotflux(2)=='fb') fradWoodf(k,i)=fWoodf(k,i)
        if(plotflux(2)=='fb') bradWoodf(k,i)=bWoodf(k,i)
      enddo
    enddo
    deallocate(Woodt)
    deallocate(Woodr)
    deallocate(Wooda)
    if(plotflux(2)=='tot') deallocate(Woodf)
    if(plotflux(2)=='fb')  deallocate(fWoodf)
    if(plotflux(2)=='fb')  deallocate(bWoodf)
    deallocate(binmaxind)
    deallocate(binmaxes)
    deallocate(fbinmax)
    deallocate(bbinmax)

if(print=='yes') print *,
if(print=='yes') print *,
if(print=='yes') print *,
if(print=='yes') print *,
if(print=='yes') print *,"radWood_rej(1)",radWood_rej(1),"    radWood_rej(2)",radWood_rej(2)
if(print=='yes') print *,
if(print=='yes') print *,"radWood refl  :",real(radWoodr(j),8)/numParts,"   radWoodr(j):",radWoodr(j)
if(print=='yes') print *,"radWood trans :",real(radWoodt(j),8)/numParts,"   radWoodt(j):",radWoodt(j)
if(print=='yes') print *,"radWood abs   :",real(radWooda(j),8)/numParts,"   radWooda(j):",radWooda(j)

  endif
  if(j==KLrnumRealz .AND. Wood=='KL') then
    allocate(KLWoodt(KLrnumRealz))
    allocate(KLWoodr(KLrnumRealz))
    allocate(KLWooda(KLrnumRealz))
    KLWood_rej(1)=Wood_rej(1)
    KLWood_rej(2)=Wood_rej(2)
    do k=1,KLrnumRealz
      KLWoodt(k)=Woodt(k)
      KLWoodr(k)=Woodr(k)
      KLWooda(k)=Wooda(k)
      do i=1,pfnumcells
        if(plotflux(2)=='tot') KLWoodf(k,i)=Woodf(k,i)
        if(plotflux(2)=='fb') fKLWoodf(k,i)=fWoodf(k,i)
        if(plotflux(2)=='fb') bKLWoodf(k,i)=bWoodf(k,i)
      enddo
    enddo
    deallocate(Woodt)
    deallocate(Woodr)
    deallocate(Wooda)
    if(plotflux(2)=='tot') deallocate(Woodf)
    if(plotflux(2)=='fb') deallocate(fWoodf)
    if(plotflux(2)=='fb') deallocate(bWoodf)
  endif


  call cpu_time(tt2)
  if(Wood=='rad') time(3)=time(3)+(tt2-tt1)
  if(Wood=='KL')  time(7)=time(7)+(tt2-tt1)

  end subroutine WoodcockMC






  subroutine WoodcockMCoutstats( numParts,radWoodt,radWoodr,&
                                 radWooda,radWood_rej,&
                                 plotflux,pltflux,pfnumcells,fluxfaces,radWoodf,&
                                 fradWoodf,bradWoodf,P )
  use genRealzvars, only: Adamscase, numRealz
  use KLvars, only: numEigs
  integer :: numParts,pfnumcells
  integer :: radWood_rej(2)
  real(8) :: P(2)
  real(8) :: fluxfaces(:),radWoodf(:,:),radWoodt(:),radWoodr(:),radWooda(:)
  real(8) :: fradWoodf(:,:),bradWoodf(:,:)
  real(8),allocatable :: Woodfave(:),Woodfvar(:),fluxinput(:)
  character(6) :: plotflux(2)
  character(7) :: pltflux(4)

  integer :: i,j
  real(8) :: radWoodr_m,radWoodr_v,radWoodt_m,radWoodt_v,radWooda_m,radWooda_v
  real(8),allocatable :: Woodfvarmat1(:),Woodfvarmat2(:)
  real(8),allocatable :: Woodfavemat1(:),Woodfavemat2(:)

  radWoodr = radWoodr/numParts
  radWoodt = radWoodt/numParts
  radWooda = radWooda/numParts

  call mean_and_var_s( radWoodr,numRealz,radWoodr_m,radWoodr_v )
  call mean_and_var_s( radWoodt,numRealz,radWoodt_m,radWoodt_v )
  call mean_and_var_s( radWooda,numRealz,radWooda_m,radWooda_v )

  633 format("radWood reflec    :   ",f8.5," +-    ",f8.5," +-  ",f10.7)
  634 format("radWood transmit  :   ",f8.5," +-    ",f8.5," +-  ",f10.7)
  635 format("radWood absorb    :   ",f8.5," +-    ",f8.5," +-  ",f10.7)

  write(*,*)
  write(*,*) "--- Transport: radWood ---"
  write(*,*) "                        mean           stdev         relerr   "
  write(*,633) radWoodr_m,sqrt(radWoodr_v),sqrt(radWoodr_v)/radWoodr_m
  write(*,634) radWoodt_m,sqrt(radWoodt_v),sqrt(radWoodt_v)/radWoodt_m
  write(*,635) radWooda_m,sqrt(radWooda_v),sqrt(radWooda_v)/radWooda_m

  636 format(" interactions     :  ",f5.1," % accept  ",f5.2," ac/p    ",f5.2," rj/p")

  write(*,636) real(radWood_rej(1),8)/(radWood_rej(1)+radWood_rej(2))*100,&
               real(radWood_rej(1),8)/numParts/numRealz,&
               real(radWood_rej(2),8)/numParts/numRealz
  write(*,*)


  open(unit=14, file="radWood.txt")    !print to data file to plot later
  640 format(f5.2,"   ",f10.4,f10.4,f10.4,f10.4,"   ",f5.2)

  write(14,*) "#Adamscase refl      refldev   trans     transdev  plot"
  write(14,640) Adamscase+0.02d0,radWoodr_m,sqrt(radWoodr_v),radWoodt_m,sqrt(radWoodt_v),Adamscase+0.08d0
  close(unit=14)
  call system("mv radWood.txt plots")


  if(pltflux(1)/='noplot') then
    if(plotflux(1)=='cell') then
      allocate(Woodfvar(pfnumcells))
      allocate(Woodfave(pfnumcells))
      allocate(fluxinput(numRealz))

      if(plotflux(2)=='tot') then
        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = radWoodf(j,i)/numParts  
          enddo
          call mean_and_var_s( fluxinput,numRealz,Woodfave(i),Woodfvar(i) )
        enddo
      endif

      if(plotflux(2)=='fb') then
        allocate(Woodfvarmat1(pfnumcells))
        allocate(Woodfvarmat2(pfnumcells))
        allocate(Woodfavemat1(pfnumcells))
        allocate(Woodfavemat2(pfnumcells))

        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = fradWoodf(j,i)/numParts/P(1)  
          enddo
          call mean_and_var_s( fluxinput,numRealz,Woodfavemat1(i),Woodfvarmat1(i) )
        enddo

        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = bradWoodf(j,i)/numParts/P(2)  
          enddo
          call mean_and_var_s( fluxinput,numRealz,Woodfavemat2(i),Woodfvarmat2(i) )
        enddo

        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = fradWoodf(j,i)/numParts + bradWoodf(j,i)/numParts  
          enddo
          call mean_and_var_s( fluxinput,numRealz,Woodfave(i),Woodfvar(i) )
        enddo
      endif


      660 format("#cell center,       ave flux,        flux dev")
      661 format(f15.7,f15.7,f15.7)
      open(unit=24, file="radWood_cellflux.txt")
      write(24,660)
      do i=1,pfnumcells
        write(24,661) (fluxfaces(i+1)+fluxfaces(i))/2.0d0,Woodfave(i),sqrt(Woodfvar(i))
      enddo
      close(unit=24)

      if(plotflux(2)=='fb') then
        680 format("#cell center,    ave flux mat1,   flux dev mat1,    ave flux mat2,   flux dev mat2")
        681 format(f15.7,f15.7,f15.7,f15.7,f15.7)
        open(unit=25, file="radWood_fbcellflux.txt")
        write(25,680)
        do i=1,pfnumcells
          write(25,681) (fluxfaces(i+1)+fluxfaces(i))/2.0d0,Woodfavemat1(i),sqrt(Woodfvarmat1(i)),&
                                                          Woodfavemat2(i),sqrt(Woodfvarmat2(i))
        enddo    
        close(unit=25)
      endif

      deallocate(Woodfvar)
      deallocate(Woodfave)
      deallocate(fluxinput)
      if(plotflux(2)=='fb') then
        deallocate(Woodfvarmat1)
        deallocate(Woodfvarmat2)
        deallocate(Woodfavemat1)
        deallocate(Woodfavemat2)
      endif

      call system("mv radWood_cellflux.txt plots")
      if(plotflux(2)=='fb') call system("mv radWood_fbcellflux.txt plots")
    endif
  endif


  end subroutine WoodcockMCoutstats 





  subroutine WoodcockKLoutstats( numParts,KLWoodt,KLWoodr,&
                                 KLWooda,KLWood_rej,&
                                 plotflux,pltflux,pfnumcells,fluxfaces,KLWoodf,&
                                 fKLWoodf,bKLWoodf,P )
  use genRealzvars, only: Adamscase
  use KLvars,       only: KLvarcalc, varmain, numEigs, KLrnumRealz
  integer :: numParts,pfnumcells
  integer :: KLWood_rej(2)
  real(8) :: P(2)
  real(8) :: fluxfaces(:),KLWoodf(:,:),KLWoodt(:),KLWoodr(:),KLWooda(:)
  real(8) :: fKLWoodf(:,:),bKLWoodf(:,:)
  real(8),allocatable :: Woodfave(:),Woodfvar(:),fluxinput(:)
  character(6) :: plotflux(2)
  character(7) :: pltflux(4)

  integer :: i,j
  real(8) :: KLWoodr_m,KLWoodr_v,KLWoodt_m,KLWoodt_v,KLWooda_m,KLWooda_v
  real(8),allocatable :: Woodfvarmat1(:),Woodfvarmat2(:)
  real(8),allocatable :: Woodfavemat1(:),Woodfavemat2(:)

  KLWoodr = KLWoodr/numParts
  KLWoodt = KLWoodt/numParts
  KLWooda = KLWooda/numParts

  call mean_and_var_s( KLWoodr,KLrnumRealz,KLWoodr_m,KLWoodr_v )
  call mean_and_var_s( KLWoodt,KLrnumRealz,KLWoodt_m,KLWoodt_v )
  call mean_and_var_s( KLWooda,KLrnumRealz,KLWooda_m,KLWooda_v )

  663 format("KLWood  reflec    :   ",f8.5," +-    ",f8.5," +-  ",f10.7)
  664 format("KLWood  transmit  :   ",f8.5," +-    ",f8.5," +-  ",f10.7)
  665 format("KLWood  absorb    :   ",f8.5," +-    ",f8.5," +-  ",f10.7)

  write(*,*)
  if(KLvarcalc=='no')  write(*,*) "--- Transport: KLWood ---"
  667 format("--- Transport: KLWood ---    Eigs=",i4,"/",i6,"    KLvar: ",f4.2)
  if(KLvarcalc=='yes') write(*,667) numEigs,size(varmain),varmain(numEigs)
  write(*,*) "                        mean           stdev         relerr   "
  write(*,663) KLWoodr_m,sqrt(KLWoodr_v),sqrt(KLWoodr_v)/KLWoodr_m
  write(*,664) KLWoodt_m,sqrt(KLWoodt_v),sqrt(KLWoodt_v)/KLWoodt_m
  write(*,665) KLWooda_m,sqrt(KLWooda_v),sqrt(KLWooda_v)/KLWooda_m

  666 format(" interactions     :  ",f5.1," % accept  ",f5.2," ac/p    ",f5.2," rj/p")

  write(*,666) real(KLWood_rej(1),8)/(KLWood_rej(1)+KLWood_rej(2))*100,&
               real(KLWood_rej(1),8)/numParts/KLrnumRealz,&
               real(KLWood_rej(2),8)/numParts/KLrnumRealz
  write(*,*)


  open(unit=14, file="KLWood.txt")    !print to data file to plot later
  640 format(f5.2,"   ",f10.4,f10.4,f10.4,f10.4,"   ",f5.2)

  write(14,*) "#Adamscase refl      refldev   trans     transdev  plot"
  write(14,640) Adamscase+.02,KLWoodr_m,sqrt(KLWoodr_v),KLWoodt_m,sqrt(KLWoodt_v),Adamscase+.08
  close(unit=14)
  call system("mv KLWood.txt plots")


  !User please be aware that fKLWoodf and bKLWoodf have no real meaning here.
  !The two materials cannot be separated because of the reconstructions.
  !The same format has been followed here as for fradWoodf and bradWoodf,
  !although only radWoodf, and not the components have any meaning.
  if(pltflux(1)/='noplot') then
    if(plotflux(1)=='cell') then
      allocate(Woodfvar(pfnumcells))
      allocate(Woodfave(pfnumcells))
      allocate(fluxinput(KLrnumRealz))

      if(plotflux(2)=='tot') then
        do i=1,pfnumcells
          do j=1,KLrnumRealz
            fluxinput(j) = KLWoodf(j,i)/numParts
          enddo
          call mean_and_var_s( fluxinput,KLrnumRealz,Woodfave(i),Woodfvar(i) )
        enddo
      endif

      if(plotflux(2)=='fb') then
        allocate(Woodfvarmat1(pfnumcells))
        allocate(Woodfvarmat2(pfnumcells))
        allocate(Woodfavemat1(pfnumcells))
        allocate(Woodfavemat2(pfnumcells))

        do i=1,pfnumcells
          do j=1,KLrnumRealz
            fluxinput(j) = fKLWoodf(j,i)/numParts/P(1)
          enddo
          call mean_and_var_s( fluxinput,KLrnumRealz,Woodfavemat1(i),Woodfvarmat1(i) )
        enddo
        do i=1,pfnumcells
          do j=1,KLrnumRealz
            fluxinput(j) = bKLWoodf(j,i)/numParts/P(2)
          enddo
          call mean_and_var_s( fluxinput,KLrnumRealz,Woodfavemat2(i),Woodfvarmat2(i) )
        enddo
        do i=1,pfnumcells
          do j=1,KLrnumRealz
            fluxinput(j) = fKLWoodf(j,i)/numParts + bKLWoodf(j,i)/numParts
          enddo
          call mean_and_var_s( fluxinput,KLrnumRealz,Woodfave(i),Woodfvar(i) )
        enddo
      endif

      670 format("#cell center,       ave flux,        flux dev")
      671 format(f15.7,f15.7,f15.7)
      open(unit=26, file="KLWood_cellflux.txt")
      write(26,670)
      do i=1,pfnumcells
        write(26,671) (fluxfaces(i+1)+fluxfaces(i))/2.0d0,Woodfave(i),sqrt(Woodfvar(i))
      enddo
      close(unit=26)

! The reason for this being commented out is at the beginning of this half of this module.
! Essentially, the results are nonsensical!
!      if(plotflux(2)=='fb') then
!        690 format("#cell center,    ave flux mat1,   flux dev mat1,    ave flux mat2,   flux dev mat2")
!        691 format(f15.7,f15.7,f15.7,f15.7,f15.7)
!        open(unit=27, file="KLWood_fbcellflux.txt")
!        write(27,690)
!        do i=1,pfnumcells
!          write(27,691) (fluxfaces(i+1)+fluxfaces(i))/2.0d0,Woodfavemat1(i),sqrt(Woodfvarmat1(i)),&
!                                                          Woodfavemat2(i),sqrt(Woodfvarmat2(i))
!        enddo    
!        close(unit=27)
!      endif

      deallocate(Woodfvar)
      deallocate(Woodfave)
      deallocate(fluxinput)
      if(plotflux(2)=='fb') then
        deallocate(Woodfvarmat1)
        deallocate(Woodfvarmat2)
        deallocate(Woodfavemat1)
        deallocate(Woodfavemat2)
      endif

      call system("mv KLWood_cellflux.txt plots")
!      if(plotflux(2)=='fb') call system("mv KLWood_fbcellflux.txt plots")!nonsensical
    endif
  endif

  end subroutine WoodcockKLoutstats





  subroutine Woodnegstats( numpnSamp,areapnSamp,distneg )
  use genRealzvars, only: numRealz
  use KLvars, only: negcnt
  integer :: numpnSamp(2)
  real(8) :: areapnSamp(4) !tot pos, tot neg, max pos, max neg
  character(3) :: distneg

  real(8) :: pos,neg

  if(distneg=='no')  print *,"--Negative smoothing stats, neg smoothing off--"
  if(distneg=='yes') print *,"--Negative smoothing stats, neg smoothing on--"
  600 format("  Neg realz   : ",f8.5,"%, ",i21," /",i21)
  601 format("  Neg samples : ",f8.5,"%, ",i21," /",i21)
  602 format("  Neg area    : ",f8.5,"%")
  603 format("  Ave neg samp: ",f11.4,"   Ave pos samp: ",f11.4)
  604 format("  Max neg samp: ",f11.4,"   Max pos samp: ",f11.4)

  write(*,600) real(negcnt,8)/real(numRealz,8),negcnt,numRealz
  pos = real(numpnSamp(1),8)
  neg = real(numpnSamp(2),8)
  write(*,601) neg/(pos+neg),numpnSamp(2),numpnSamp(1)+numpnSamp(2)
  write(*,602) -areapnSamp(2)/(areapnSamp(1)-areapnSamp(2))
  write(*,603) areapnSamp(2)/numpnSamp(2),areapnSamp(1)/numpnSamp(1)
  write(*,604) areapnSamp(4),areapnSamp(3)
  end subroutine





  subroutine transplot( Adamscase )
  real(8) :: Adamscase

  open(unit=17, file="tempAdams")                                  !create script
  650 format("grep '",f3.1,"         ' 'texts/transAdams_plot.txt' > plots/transAdams.txt")
  write(17,650) Adamscase
  close(unit=17)
  call system("chmod u+x tempAdams")
  call system("./tempAdams")                                       !use script
  call system("rm tempAdams")                                      !discard script

  call system("gnuplot plots/transportgnus/transport.gnu")         !create plot
  call system("ps2pdf transport.ps transport.pdf")                 !other formats
  call system("ps2eps transport.ps")
  call system("mv transport.ps transport.pdf transport.eps plots") !archive

  end subroutine transplot



  subroutine KLWood_binmaxes( j,lamc,binmaxind,binmaxes,nbin)
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave
  integer :: j,nbin
  real(8) :: lamc
  real(8) :: binmaxind(:),binmaxes(:)

  integer :: i,k
  integer :: numinnersteps = 7
  integer :: numrefine     = 5
  real(8) :: safetyfactor  = 1.1d0
  real(8) :: innerstep,maxsig,maxpos,xpos,xsig,xpos1,xsig1,xpos2,xsig2

  do i=1,nbin
    innerstep = (binmaxind(2)-binmaxind(1)) / (numinnersteps-1)
    maxsig=0.0d0             !find initial max val
    maxpos=0.0d0
    do k=1,numinnersteps
      xpos=binmaxind(i)+(k-1)*innerstep
      xsig= KLrxi_point(j,lamc,xpos)
      if(xsig>maxsig) then
        maxsig=xsig
        maxpos=xpos
      endif
    enddo
    do k=1,numrefine     !refine
      innerstep=innerstep/2
      xpos1=maxpos-innerstep
      xsig1= KLrxi_point(j,lamc,xpos1)
      xpos2=maxpos+innerstep
      xsig2= KLrxi_point(j,lamc,xpos2)
      if(xsig1>maxsig .AND. xsig1>xsig2) then
        maxsig=xsig1
        maxpos=xpos1
      elseif(xsig2>maxsig .AND. xsig2>xsig1) then
        maxsig=xsig2
        maxpos=xpos2
      else
        maxpos=xpos
      endif
    enddo
    binmaxes(i)=maxsig*safetyfactor   !assign maxsig

  enddo

  end subroutine KLWood_binmaxes




  function ceilsigfunc(binmaxind,binmax,position,nbin)
  integer :: nbin
  real(8) :: binmaxind(:),binmax(:),position,ceilsigfunc

  integer :: i

  do i=1,nbin
    if( binmaxind(i)<=position .AND. binmaxind(i+1)>position ) then
      ceilsigfunc = binmax(i)
      exit
    endif
  enddo

  end function ceilsigfunc


  function radWood_actsig(position,matType,matLength,sig)
  integer :: matType(:)
  real(8) :: position,matLength(:),sig(2),radWood_actsig

  integer :: i

  i=1
  do
    if(matLength(i)<=position .AND. matLength(i+1)>position) then
      radWood_actsig=sig(matType(i))
      exit
    endif 
    i=i+1
  enddo

  end function radWood_actsig


  function radWood_actscatrat(position,matType,matLength,scatrat)
  integer :: matType(:)
  real(8) :: position,matLength(:),scatrat(2),radWood_actscatrat

  integer :: i

  i=1
  do
    if(matLength(i)<=position .AND. matLength(i+1)>position) then
      radWood_actscatrat=scatrat(matType(i))
      exit
    endif
    i=i+1
  enddo

  end function radWood_actscatrat


  subroutine radWood_binmaxes(matLength,matType,numArrSz,binmaxind,binmaxes,nbin,sig)
  !subroutine starts to set up ceiling for WoodcockMC by mapping highest point in each bin
  integer :: numArrSz,nbin,matType(:)
  real(8) :: matLength(:),binmaxind(:),binmaxes(:),sig(2)

  integer :: i,k
  real(8) :: smallersig,largersig

  smallersig=minval(sig)
  largersig =maxval(sig)

  do i=1,nbin
    do k=1,numArrSz

      !contains both
      if(matLength(k)>binmaxind(i) .AND. matLength(k)<binmaxind(i+1)) then
        binmaxes(i)=largersig
        exit
      endif

      !contains only one
      if(matLength(k)<=binmaxind(i) .AND. matLength(k+1)>=binmaxind(i+1)) then
        binmaxes(i)=sig(matType(k))
        exit
      endif

    enddo
  enddo

  end subroutine





end module

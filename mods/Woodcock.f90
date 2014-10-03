module Woodcock
  use mcnp_random
  use timeman
  use radtransMC
  use utilities
  use KLmeanadjust
  implicit none

CONTAINS
  ! print statements in this module use 600-699

  subroutine WoodcockMC( j )
  use timevars, only: time
  use genRealzvars, only: sig, scatrat, lam, s, numRealz, nummatSegs, lamc, &
                          matType, matLength
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave, KLrnumRealz
  use MCvars, only: numParts, pfnumcells, rodOrplanar, sourceType, plotflux, &
                    pltflux, Woodt, Woodr, radWoodr, KLWoodr, radWoodt, KLWoodt, &
                    Wooda, radWooda, KLWooda, Wood_rej, radWood_rej, KLWood_rej, &
                    numpnSamp, areapnSamp, allowneg, distneg, Wood, fluxfaces, &
                    fWoodf, bWoodf, fradWoodf, bradWoodf, fKLWoodf, bKLWoodf, &
                    radWoodf, KLWoodf, Woodf, binmaxind, binmaxes, fbinmax, bbinmax, &
                    nceilbin
  integer :: j

  integer :: i,o,k
  real(8) :: disthold
  real(8) :: binlength,position,mu,ceilsig,woodrat,actsig,tscatrat
  real(8) :: db,dc,tt1,tt2

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
  if(Wood=='rad') nceilbin = ceiling(s/lamc)
  if(Wood=='KL')  nceilbin = numEigs

  allocate(binmaxind(nceilbin+1))
  allocate(binmaxes(nceilbin))
  allocate(fbinmax(nceilbin))
  allocate(bbinmax(nceilbin))
  binmaxind=0.0d0
  binmaxes=0.0d0
  fbinmax=0.0d0
  bbinmax=0.0d0

  binlength=s/nceilbin
  binmaxind(1)=0.0d0
  do i=2,nceilbin+1
    binmaxind(i)=binmaxind(i-1)+binlength
    if(print=='yes') print *,"i:",i,"binmaxind(i):",binmaxind(i),"   nceilbin",nceilbin
  enddo

  if(Wood=='rad') call radWood_binmaxes( nummatSegs )
  if(Wood=='KL')  call KLWood_binmaxes( j )

  !create forward/backward motion max vectors
  bbinmax(1)=binmaxes(1)
  fbinmax(nceilbin)=binmaxes(nceilbin)
  do i=2,nceilbin
    bbinmax(i) = merge(bbinmax(i-1),binmaxes(i),bbinmax(i-1)>binmaxes(i))
  enddo
  do i=nceilbin-1,1,-1
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
      
      ceilsig=merge(ceilsigfunc(position,fbinmax),& !sel max sig
             ceilsigfunc(position,bbinmax),mu>=0)

      dc = -log(1-rang())/ceilsig                   !calc dc

      if(print=='yes') print *,"  position   :",position,"    mu     :",mu,&
                               "   db:",db,"   ceilsig:",ceilsig,"   dc:",dc


      if( dc>db ) Wood_rej(1)=Wood_rej(1)+1   !accept path
      if( dc>db .AND. mu>=0 ) then        !tally trans or refl
        Woodt(j) = Woodt(j)+1 !transmit
        if(plotflux(2)=='tot') call adv_pos_col_flux(position,s,Woodf,j,mu)
        if(plotflux(2)=='fb') call col_fbflux(position,s,fWoodf,bWoodf,j,mu)

        if(print=='yes') print *,"                      tally transmit"
        exit
      endif
      if( dc>db .AND. mu<0 ) then
        Woodr(j) = Woodr(j) + 1 !reflect
        if(plotflux(2)=='tot') call adv_pos_col_flux(position,0.0d0,Woodf,j,mu)
        if(plotflux(2)=='fb') call col_fbflux(position,0.0d0,fWoodf,bWoodf,j,mu)
        if(print=='yes') print *,"                      tally reflect"
        exit
      endif

      if(plotflux(2)=='tot') call adv_pos_col_flux(position,position+dc*mu,Woodf,j,mu)
      if(plotflux(2)=='fb') call col_fbflux(position,position+dc*mu,fWoodf,bWoodf,j,mu)
      if(Wood=='rad') woodrat= radWood_actsig(position,sig)/ceilsig
      if(Wood=='KL')  woodrat= KLrxi_point(j,position)/ceilsig
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

      if(Wood=='rad') tscatrat=radWood_actscatrat(position,scatrat)
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









  subroutine Woodnegstats
  use genRealzvars, only: numRealz
  use KLvars, only: negcnt
  use MCvars, only: numpnSamp, areapnSamp, distneg, KLWood, allowneg

  real(8) :: pos,neg

  open(unit=100,file="Woodnegstats.out")

  if(KLWood=='yes' .and. allowneg=='yes') then
    if(distneg=='no')  write(100,*) "--Negative smoothing stats, neg smoothing off--"
    if(distneg=='yes') write(100,*) "--Negative smoothing stats, neg smoothing on--"
    600 format("  Neg realz   : ",f8.5,"%, ",i21," /",i21)
    601 format("  Neg samples : ",f8.5,"%, ",i21," /",i21)
    602 format("  Neg area    : ",f8.5,"%")
    603 format("  Ave neg samp: ",f11.4,"   Ave pos samp: ",f11.4)
    604 format("  Max neg samp: ",f11.4,"   Max pos samp: ",f11.4)

    write(100,600) real(negcnt,8)/real(numRealz,8),negcnt,numRealz
    pos = real(numpnSamp(1),8)
    neg = real(numpnSamp(2),8)
    write(100,601) neg/(pos+neg),numpnSamp(2),numpnSamp(1)+numpnSamp(2)
    write(100,602) -areapnSamp(2)/(areapnSamp(1)-areapnSamp(2))
    write(100,603) areapnSamp(2)/numpnSamp(2),areapnSamp(1)/numpnSamp(1)
    write(100,604) areapnSamp(4),areapnSamp(3)
    write(100,*)
  else
    write(100,*)
  endif

  close(unit=100)
  call system("mv Woodnegstats.out texts")
  end subroutine Woodnegstats







  subroutine KLWood_binmaxes( j )
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave
  use MCvars, only: binmaxind, binmaxes, nceilbin
  integer :: j

  integer :: i,k
  integer :: numinnersteps = 7
  integer :: numrefine     = 5
  real(8) :: safetyfactor  = 1.1d0
  real(8) :: innerstep,maxsig,maxpos,xpos,xsig,xpos1,xsig1,xpos2,xsig2

  do i=1,nceilbin
    innerstep = (binmaxind(2)-binmaxind(1)) / (numinnersteps-1)
    maxsig=0.0d0             !find initial max val
    maxpos=0.0d0
    do k=1,numinnersteps
      xpos=binmaxind(i)+(k-1)*innerstep
      xsig= KLrxi_point(j,xpos)
      if(xsig>maxsig) then
        maxsig=xsig
        maxpos=xpos
      endif
    enddo
    do k=1,numrefine     !refine
      innerstep=innerstep/2
      xpos1=maxpos-innerstep
      xsig1= KLrxi_point(j,xpos1)
      xpos2=maxpos+innerstep
      xsig2= KLrxi_point(j,xpos2)
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




  function ceilsigfunc(position,binmax)
  use MCvars, only: nceilbin, binmaxind
  real(8) :: position,ceilsigfunc,binmax(:)

  integer :: i

  do i=1,nceilbin
    if( binmaxind(i)<=position .AND. binmaxind(i+1)>position ) then
      ceilsigfunc = binmax(i)
      exit
    endif
  enddo

  end function ceilsigfunc


  function radWood_actsig(position,sig)
  use genRealzvars, only: matType, matLength
  real(8) :: position,sig(2),radWood_actsig

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


  function radWood_actscatrat(position,scatrat)
  use genRealzvars, only: matType, matLength
  real(8) :: position,scatrat(2),radWood_actscatrat

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


  subroutine radWood_binmaxes( numArrSz )
  !subroutine starts to set up ceiling for WoodcockMC by mapping highest point in each bin
  use genRealzvars, only: matType, matLength, sig
  use MCvars, only: nceilbin, binmaxind, binmaxes
  integer :: numArrSz

  integer :: i,k
  real(8) :: smallersig,largersig

  smallersig=minval(sig)
  largersig =maxval(sig)

  do i=1,nceilbin
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

  end subroutine radWood_binmaxes





end module

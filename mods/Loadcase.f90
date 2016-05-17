module Loadcase
  implicit none

CONTAINS
  ! print statements in this module use # 100-199

  subroutine read_test_inputstoc
  use rngvars, only: rngseed
  use genRealzvars,         only: Adamscase, sig, scatrat, lam, slen, numRealz, pltgenrealznumof, &
                                  pltgenrealz, pltgenrealzwhich, GBaves1, GBavea1, GBaves2, GBavea2, &
                                  GBvars1, GBvara1, GBvars2, GBvara2, &
                                  GBlamcs1,GBlamca1,GBlamcs2,GBlamca2, chgeomtype
  use KLvars,               only: pltEigfwhich, pltxiBinswhich, numEigss1, numEigsa1, numEigss2, numEigsa2,&
                                  pltCowhich, pltxiBinsnumof, pltEigfnumof, pltConumof, binNumof,&
                                  numSlice, levsrefEig, Corrnumpoints, binSmallBound, &
                                  binLargeBound, pltxiBins, pltxiBinsgauss, pltEigf, pltCo, &
                                  Corropts, KLrnumpoints, pltKLrealz, pltKLrealznumof, pltKLrealzwhich, &
                                  flmeanadjust, meanadjust_tol, chGBcase, corrind, &
                                  Gaussrandtype, numrefinesameiter, corrinds1, corrinda1, corrinds2, &
                                  corrinda2, chGausstype, lamctypes1, lamctypea1, lamctypes2, lamctypea2, &
                                  numLNxspts, numLNxsbins, chLNxschecktype, chLNxsplottype, &
                                  fls1, fla1, fls2, fla2, numNystroms1, numNystroma1, numNystroms2, numNystroma2, &
                                  cheftypes1, cheftypea1, cheftypes2, cheftypea2
  use MCvars,               only: trprofile_binnum, binplot, numParts, trannprt, rodOrplanar, sourceType, &
                                  pltflux, flnegxs, LPamnumParts, fluxnumcells, pltmatflux, mindatapts, &
                                  pltfluxtype, flCR_MCSC, chTrantype, reflrelSEMtol, tranrelSEMtol, maxnumParts
  use UQvars,               only: Qs, numUQdims, chUQtype, PCEorder, flPCErefl, flPCEtran, PCEcells, numPCEcells, &
                                  flPCEsampscorrelated, numPCEQoIsamps
  use rngvars,              only: rngstride
  character(7) :: pltallopt                         !Plot all same opt

  character(4)  :: setflags(3)
  character(20) :: dumchar !use this to "skip" a line
  integer         :: i, ic, Qtempsize, valsused, maxnumeigs
  integer, allocatable :: Qtemp(:),tcorrind(:),Qtemparray(:,:)
  logical         :: flstopstatus = .false., flsleep = .false.

  open(unit=2,file="inputstoc.txt")

  read(2,*) rngseed

  !--- Biggest Problem Parameters ---!
  read(2,*) dumchar
  read(2,*) chgeomtype
  read(2,*) chTrantype
  read(2,*) chUQtype
  read(2,*) numRealz,trannprt
  read(2,*) numParts,maxnumParts,reflrelSEMtol,tranrelSEMtol,mindatapts
  read(2,*) slen

  !--- Geometry - Gauss or Gauss-based type problem ---!
  read(2,*) dumchar
  read(2,*) chGausstype,chGBcase
  read(2,*) fls1,GBaves1,GBvars1,GBlamcs1,numEigss1,lamctypes1,corrinds1,numNystroms1,cheftypes1
  read(2,*) fla1,GBavea1,GBvara1,GBlamca1,numEigsa1,lamctypea1,corrinda1,numNystroma1,cheftypea1
  read(2,*) fls2,GBaves2,GBvars2,GBlamcs2,numEigss2,lamctypes2,corrinds2,numNystroms2,cheftypes2
  read(2,*) fla2,GBavea2,GBvara2,GBlamca2,numEigsa2,lamctypea2,corrinda2,numNystroma2,cheftypea2
  if(chUQtype=="SC" .or. chUQtype=="PCE") then !test KL orders when using cubature
    if(fls1 .and. fla1 .and. abs(corrinds1)==abs(corrinda1) .and. (numEigss1/=numEigsa1)) flstopstatus = .true.
    if(fls1 .and. fls2 .and. abs(corrinds1)==abs(corrinds2) .and. (numEigss1/=numEigss2)) flstopstatus = .true.
    if(fls1 .and. fla2 .and. abs(corrinds1)==abs(corrinda2) .and. (numEigss1/=numEigsa2)) flstopstatus = .true.
    if(fla1 .and. fls2 .and. abs(corrinda1)==abs(corrinds2) .and. (numEigsa1/=numEigss2)) flstopstatus = .true.
    if(fla1 .and. fla2 .and. abs(corrinda1)==abs(corrinda2) .and. (numEigsa1/=numEigsa2)) flstopstatus = .true.
    if(fls2 .and. fla2 .and. abs(corrinds2)==abs(corrinda2) .and. (numEigss2/=numEigsa2)) flstopstatus = .true.
    if(flstopstatus) print *,"--For corr or anticorr KL expansions with SC or PCE, order must be the same"
    if(flstopstatus) stop
  endif

  if(fls1) then              !Test for approriate correlation indices (first==1, any other
    allocate(corrind(1))     !number preceded by positve counterpart or preceding integer.  
    corrind = corrinds1
  endif
  if(fla1) then
    if(allocated(corrind)) then
      call move_alloc(corrind,tcorrind)
      allocate(corrind(size(tcorrind)+1))
      corrind(:size(tcorrind)) = tcorrind
      deallocate(tcorrind)
    else
      allocate(corrind(1))
    endif
    corrind(size(corrind)) = corrinda1
  endif
  if(fls2) then
    if(allocated(corrind)) then
      call move_alloc(corrind,tcorrind)
      allocate(corrind(size(tcorrind)+1))
      corrind(:size(tcorrind)) = tcorrind
      deallocate(tcorrind)
    else
      allocate(corrind(1))
    endif
    corrind(size(corrind)) = corrinds2
  endif
  if(fla2) then
    if(allocated(corrind)) then
      call move_alloc(corrind,tcorrind)
      allocate(corrind(size(tcorrind)+1))
      corrind(:size(tcorrind)) = tcorrind
      deallocate(tcorrind)
    else
      allocate(corrind(1))
    endif
    corrind(size(corrind)) = corrinda2
  endif

  if(corrind(1)/=1) then
    print *,"--First used correlation indice must be '1'"
    flstopstatus = .true.
  endif
  do i=1,size(corrind)
    if(i>1) then
      if(corrind(i)<0) then
        if(.not. any(corrind(:i)==abs(corrind(i)))) then
          print *,"--Negative correlation indices must be preceded by positive counterpart"
          flstopstatus = .true.
        endif
      endif
      if(corrind(i)>1) then
        if(.not. any(corrind(:i)==corrind(i)-1)) then
          print *,"--Positive correlation indices must be preceded by previous integer"
          flstopstatus = .true.
        endif
      endif
    endif
  enddo

  if(chUQtype=="SC" .or. chUQtype=="PCE") then
    numUQdims = 0 !cycle through correlations and count number of KL variables
    do ic=1,4
      if(fls1 .and. ic==abs(corrinds1)) then
        numUQdims = numUQdims + numEigss1
        cycle
      endif
      if(fla1 .and. ic==abs(corrinda1)) then
        numUQdims = numUQdims + numEigsa1
        cycle
      endif
      if(fls2 .and. ic==abs(corrinds2)) then
        numUQdims = numUQdims + numEigss2
        cycle
      endif
      if(fla2 .and. ic==abs(corrinda2)) then
        numUQdims = numUQdims + numEigsa2
        cycle
      endif
    enddo
  else
    numUQdims = 1
  endif
  allocate(Qs(numUQdims)) !allocate and read quadrature orders according to correlations
  read(2,*) (Qs(i),i=1,numUQdims)
  read(2,*) PCEorder

  !--- Geometry - 'Markov' type problem ---!
  read(2,*) dumchar
  read(2,*) Adamscase
  read(2,*) sig(1),sig(2)
  read(2,*) scatrat(1),scatrat(2)
  read(2,*) lam(1),lam(2)

  !--- Other KL Options ---!
  read(2,*) dumchar
  read(2,*) binNumof
  read(2,*) KLrnumpoints
  read(2,*) levsrefEig
  read(2,*) numrefinesameiter
  read(2,*) binSmallBound,binLargeBound
  read(2,*) numSlice
  read(2,*) Gaussrandtype

  !--- Other PCE Options ---!
  read(2,*) dumchar
  read(2,*) numPCEQoIsamps
  read(2,*) flPCEsampscorrelated
  read(2,*) flPCErefl,flPCEtran
  read(2,*) numPCEcells
  if(numPCEcells>=1) then
    allocate(PCEcells(numPCEcells))
    read(2,*) (PCEcells(i),i=1,numPCEcells)
  endif
  if(numPCEcells<1) read(2,*) dumchar !0 means no PCE cells, -1 is all, handle later

  !--- Other MCtrans Options ---!
  read(2,*) dumchar 
  read(2,*) LPamnumParts
  read(2,*) setflags(1)
  if(setflags(1)=='yes') flCR_MCSC    =.true.
  read(2,*) rodOrplanar
  read(2,*) sourceType
  read(2,*) setflags(1)
  if(setflags(1)=='yes') flnegxs  =.true.
  read(2,*) setflags(1),meanadjust_tol
  if(setflags(1)=='yes') flmeanadjust=.true.


  read(2,*) dumchar    !All Plot Same Way Option
  read(2,*) pltallopt

  read(2,*) dumchar    !Plotting genRealz realz
  read(2,*) pltgenrealz(1),pltgenrealz(2),pltgenrealz(3),pltgenrealz(4)
  read(2,*) pltgenrealznumof
  allocate(pltgenrealzwhich(pltgenrealznumof))
  do i=1,pltgenrealznumof
    read(2,*) pltgenrealzwhich(i)
  enddo


  read(2,*) dumchar    !Plotting Eigenfunction
  read(2,*) pltEigf(1),pltEigf(2),pltEigf(3),pltEigf(4)
  read(2,*) pltEigfnumof
  allocate(pltEigfwhich(pltEigfnumof))
  do i=1,pltEigfnumof
    read(2,*) pltEigfwhich(i)
  enddo

  read(2,*) dumchar    !Plotting xiBins
  read(2,*) pltxiBins(1),pltxiBins(2),pltxiBins(3),pltxiBins(4)
  read(2,*) pltxiBinsgauss
  read(2,*) pltxiBinsnumof
  allocate(pltxiBinswhich(2,pltxiBinsnumof))
  do i=1,pltxiBinsnumof
    read(2,*) pltxiBinswhich(1,i),pltxiBinswhich(2,i)
  enddo

  read(2,*) dumchar    !Plotting KLrealz
  read(2,*) pltKLrealz(1),pltKLrealz(2),pltKLrealz(3),pltKLrealz(4)
  read(2,*) pltKLrealznumof
  allocate(pltKLrealzwhich(3,pltKLrealznumof))
  do i=1,pltKLrealznumof
    read(2,*) pltKLrealzwhich(1,i),pltKLrealzwhich(2,i),pltKLrealzwhich(3,i)
  enddo

  read(2,*) dumchar    !Plotting local KL average and variance
  read(2,*) chLNxsplottype,chLNxschecktype
  read(2,*) numLNxspts
  read(2,*) numLNxsbins

  read(2,*) dumchar    !Plotting Variace (Co)
  read(2,*) pltCo(1),pltCo(2),pltCo(3),pltCo(4)
  read(2,*) pltConumof
  allocate(pltCowhich(pltConumof))
  do i=1,pltConumof
    read(2,*) pltCowhich(i)
  enddo

  read(2,*) dumchar    !Plotting Correlation contours
  read(2,*) Corropts(1),Corropts(2)
  read(2,*) Corrnumpoints


  read(2,*) dumchar    !Leakage pdf
  read(2,*) binplot
  read(2,*) trprofile_binnum

  read(2,*) dumchar    !Plotting flux
  read(2,*) pltflux(1),pltflux(2),pltflux(3),pltflux(4)
  read(2,*) pltmatflux
  read(2,*) pltfluxtype
  read(2,*) fluxnumcells


  if( pltallopt .ne. 'default' ) then
    chLNxsplottype = pltallopt
    binplot        = pltallopt
    pltEigf(1)     = pltallopt
    Corropts(1)    = pltallopt
    pltxiBins(1)   = pltallopt
    pltKLrealz(1)  = pltallopt
    pltgenrealz(1) = pltallopt
    pltCo(1)       = pltallopt
    pltflux(1)     = pltallopt
    pltmatflux     = pltallopt
  endif


  if(chgeomtype=='binary' .and. Adamscase/=0) call Acase_load !need to load these to test
  if(chgeomtype=='contin' .and. .not.chGBcase=='none') call GBcase_load


  !begin tests of valid input
  print *,"  "

  !Test PCEcells options, allocate PCEcells to all flux cells if necessary
  if(chUQtype=='PCE') then
    if(numPCEcells/=0 .and. .not.(pltflux(1)=='plot' .or. pltflux(1)=='preview')) then
      print *,"--User choosing PCE over cells without plotting flux"
      flstopstatus = .true.
    endif
    if(numPCEcells==-1) then
      allocate(PCEcells(fluxnumcells))
      numPCEcells = fluxnumcells
    elseif(numPCecells<-1) then
      print *,"--User must enter an integer in [-1,infty) for number of PCE cells"
      flstopstatus = .true.
    endif
    if(numPCEcells > fluxnumcells) then
      print *,"--User attemping to apply PCE to more cells than plotting flux for"
      flstopstatus = .true.
    endif
    if(numPCEcells>0) then
      do i=1,numPCEcells
        if(PCEcells(i)<1 .or. PCEcells(i)>fluxnumcells) then
          print *,"--User attemping to apply PCE to cell number which is negative or greater than number of flux cells"
          flstopstatus = .true.
        endif
      enddo
    endif
  endif

  !Test Nystrom options
  if(fls1 .and. lamctypes1=='numeric' .and. numEigss1>numNystroms1) then
    print *,"--User specified coarser Nymstrom discretization than number of eigenvalues, s1"
    flstopstatus = .true.
  endif
  if(fla1 .and. lamctypea1=='numeric' .and. numEigsa1>numNystroma1) then
    print *,"--User specified coarser Nymstrom discretization than number of eigenvalues, a1"
    flstopstatus = .true.
  endif
  if(fls2 .and. lamctypes2=='numeric' .and. numEigss2>numNystroms2) then
    print *,"--User specified coarser Nymstrom discretization than number of eigenvalues, s2"
    flstopstatus = .true.
  endif
  if(fla2 .and. lamctypea2=='numeric' .and. numEigsa2>numNystroma2) then
    print *,"--User specified coarser Nymstrom discretization than number of eigenvalues, a2"
    flstopstatus = .true.
  endif

  !Tests for problem type
  if(.not.chgeomtype=='contin' .and. .not.chgeomtype=='binary') then
    print *,"--User specified illegal geomtype.  Options: 'contin' or 'binary'"
    flstopstatus = .true.
  endif
  if(chgeomtype=='contin' .and. .not.(chTrantype=='GaussKL' .or. chTrantype=='None')) then
    print *,"--User specified contin geom but not GaussKL transport, set to GaussKL"
    chTrantype = 'GaussKL'
    flsleep = .true.
  endif
  if(chgeomtype=='binary' .and. .not.(chTrantype=='radMC' .or. chTrantype=='radWood' .or.&
     chTrantype=='KLWood'.or. chTrantype=='LPMC'  .or. chTrantype=='atmixMC' .or. chTrantype=='None')) then
    print *,"--User attempting to use invalid transport type for binary material geometry"
    flstopstatus = .true.
  endif
  if(chgeomtype=='contin' .and. .not.(chGBcase=='none' .or. chGBcase=='f1' .or. chGBcase=='f2' .or. &
                                                    chGBcase=='ANS16l' .or. chGBcase=='ANS16s' .or. &
                                      chGBcase=='A1.1' .or.chGBcase=='A1.2'.or.chGBcase=='A1.3'.or. &
                                      chGBcase=='A2.1' .or.chGBcase=='A2.2'.or.chGBcase=='A2.3'.or. &
                                      chGBcase=='A3.1' .or.chGBcase=='A3.2'.or.chGBcase=='A3.3'.or. &
                                      chGBcase=='A4.1' .or.chGBcase=='A4.2'.or.chGBcase=='A4.3'.or. &
                                      chGBcase=='A5.1' .or.chGBcase=='A5.2'.or.chGBcase=='A5.3'.or. &
                                      chGBcase=='A6.1' .or.chGBcase=='A6.2'.or.chGBcase=='A6.3'.or. &
                                      chGBcase=='A7.1' .or.chGBcase=='A7.2'.or.chGBcase=='A7.3'.or. &
                                      chGBcase=='A8.1' .or.chGBcase=='A8.2'.or.chGBcase=='A8.3'.or. &
                                      chGBcase=='A9.1' .or.chGBcase=='A9.2'.or.chGBcase=='A9.3'.or. &
                                      chGBcase=='A1.1two' .or.chGBcase=='A1.2two'.or.chGBcase=='A1.3two'.or. &
                                      chGBcase=='A2.1two' .or.chGBcase=='A2.2two'.or.chGBcase=='A2.3two'.or. &
                                      chGBcase=='A3.1two' .or.chGBcase=='A3.2two'.or.chGBcase=='A3.3two'.or. &
                                      chGBcase=='A4.1two' .or.chGBcase=='A4.2two'.or.chGBcase=='A4.3two'.or. &
                                      chGBcase=='A5.1two' .or.chGBcase=='A5.2two'.or.chGBcase=='A5.3two'.or. &
                                      chGBcase=='A6.1two' .or.chGBcase=='A6.2two'.or.chGBcase=='A6.3two'.or. &
                                      chGBcase=='A7.1two' .or.chGBcase=='A7.2two'.or.chGBcase=='A7.3two'.or. &
                                      chGBcase=='A8.1two' .or.chGBcase=='A8.2two'.or.chGBcase=='A8.3two'.or. &
                                      chGBcase=='A9.1two' .or.chGBcase=='A9.2two'.or.chGBcase=='A9.3two'.or. &
                                      chGBcase=='A1.1ind' .or.chGBcase=='A1.2ind'.or.chGBcase=='A1.3ind'.or. &
                                      chGBcase=='A2.1ind' .or.chGBcase=='A2.2ind'.or.chGBcase=='A2.3ind'.or. &
                                      chGBcase=='A3.1ind' .or.chGBcase=='A3.2ind'.or.chGBcase=='A3.3ind'.or. &
                                      chGBcase=='A4.1ind' .or.chGBcase=='A4.2ind'.or.chGBcase=='A4.3ind'.or. &
                                      chGBcase=='A5.1ind' .or.chGBcase=='A5.2ind'.or.chGBcase=='A5.3ind'.or. &
                                      chGBcase=='A6.1ind' .or.chGBcase=='A6.2ind'.or.chGBcase=='A6.3ind'.or. &
                                      chGBcase=='A7.1ind' .or.chGBcase=='A7.2ind'.or.chGBcase=='A7.3ind'.or. &
                                      chGBcase=='A8.1ind' .or.chGBcase=='A8.2ind'.or.chGBcase=='A8.3ind'.or. &
                                      chGBcase=='A9.1ind' .or.chGBcase=='A9.2ind'.or.chGBcase=='A9.3ind'     ) ) then
    print *,"--User giving non-valid option for special GB case"
    flstopstatus = .true.
  endif
  if(chgeomtype=='binary' .and. .not.chUQtype=='MC') then
    print *,"--User attempting to use non-MC UQ method with binary geometry"
    flstopstatus = .true.
  endif
  if(rngstride < numRealz) then
    print *,"--User using too many realizations--reusing random numbers, may cause correlation"
    flsleep = .true.
  endif
  if((chTrantype=='radMC'  .or. chTrantype=='radWood' .or. &
      chTrantype=='KLWood' .or.  chTrantype=='GaussKL') .and. maxnumParts<numParts) then
    print *,"--maxnumParts smaller than, so set equal to, numParts"
    maxnumParts = numParts
  endif
  if((chTrantype=='LPMC' .or. chTrantype=='atmixMC') .and. maxnumParts<LPamnumParts) then
    print *,"--maxnumParts smaller than, so set equal to, LPamnumParts"
    maxnumParts = LPamnumParts
  endif


  do i=1,pltEigfnumof    !Test Eigenfunction plotting order of Eigs
    if( pltEigfwhich(i)>numEigss1 .AND. pltEigf(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot Eigenfunction of higher order than Eigenvalues calculated"
      flstopstatus = .true.
    endif
  enddo

  do i=1,pltxiBinsnumof  !Test xiBins plotting order of Eigs and num of bins
    if( pltxiBinswhich(1,i)>numEigss1 .AND. pltxiBins(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot xiBins for Eigenvalue of higher order than calculated"
      flstopstatus = .true.
    endif
    if( pltxiBinswhich(2,i)>numRealz .AND. pltxiBins(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot xiBins for more realizations than generated"
      flstopstatus = .true.
    endif
  enddo

  do i=1,pltKLrealznumof
    if( pltKLrealzwhich(1,i)>numRealz .AND. pltKLrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot more reconstructed realz than reconstructed"
      flstopstatus = .true.
    endif
    if( (pltKLrealzwhich(2,i)>numEigss1 .or. pltKLrealzwhich(3,i)>numEigsa1) .and. pltKLrealz(1) .ne. 'noplot' ) then
      print *,"--User attempting to plot reconstructed realz using more than calced num of Eigs"
      flstopstatus = .true.
    endif
  enddo
  if( Gaussrandtype .ne. 'BM' .and. Gaussrandtype .ne. 'inv' ) then
    print *,"--User should enter 'BM' or 'inv' for Gaussian sampling type"
    flstopstatus = .true.
  endif
  !Test for meanadjust
  if(.not.flnegxs .and. flmeanadjust) then
    print *,"--User chose meanadjust and disallowed realz w/ neg, meanadjust turned off"
    flmeanadjust = .false.
    flsleep = .true.
  endif

                              !Test plotting flux
  if( pltflux(1)/='noplot' .and. .not.chTrantype=='radMC' .and. .not.chTrantype=='radWood' &
      .and. .not.chTrantype=='KLWood' .and. .not.chTrantype=='LPMC' .and. &
      .not.chTrantype=='atmixMC' .and. .not.chTrantype=='GaussKL') then
    print *,"--User attempting to plot flux when no transport calculations are made"
    flstopstatus = .true.
  endif
  if( sourceType/='leftbeam' .and. sourceType/='leftiso' .and. sourceType/='intern' ) then
    print *,"--User attempting to run invalid source type.  Options: 'leftbeam','leftiso','intern'"
    flstopstatus = .true.
  endif
  if( pltfluxtype/='track' .and. pltfluxtype/='point' ) then
    print *,"--User attempting to plot flux with invalid scheme.  Please enter 'track' or 'point'"
    flstopstatus = .true.
  endif


  do i=1,pltgenrealznumof    !Test genRealz plotting over selected realz
    if( pltgenrealzwhich(i)>numRealz .AND. pltgenrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot realizations that are not created"
      flstopstatus = .true.
    endif
  enddo

  do i=1,pltConumof          !Test plotCo for Eig choice and CoEffExp or CoEffAct option
    if( pltCowhich(i)>numEigss1 .and. pltCo(1) .ne. 'noplot' ) then
      print *,"--User attempting to plot CoEff values for eigenvalues not calculated"
      flstopstatus = .true.
    endif
  enddo

  !Test plotting valid options
  if(.not.chLNxsplottype=='noplot' .and. .not.chLNxsplottype=='plot' .and. .not.chLNxsplottype=='preview') then
    print *,"--Use given invalid option for chLNxsplottype"
    flstopstatus = .true.
  endif
  if(.not.binplot=='noplot' .and. .not.binplot=='plot' .and. .not.binplot=='preview') then
    print *,"--Use given invalid option for binplot"
    flstopstatus = .true.
  endif
  if(.not.pltEigf(1)=='noplot' .and. .not.pltEigf(1)=='plot' .and. .not.pltEigf(1)=='preview') then
    print *,"--Use given invalid option for pltEigf(1)"
    flstopstatus = .true.
  endif
  if(.not.Corropts(1)=='noplot' .and. .not.Corropts(1)=='plot' .and. .not.Corropts(1)=='preview') then
    print *,"--Use given invalid option for Corropts(1)"
    flstopstatus = .true.
  endif
  if(.not.pltxiBins(1)=='noplot' .and. .not.pltxiBins(1)=='plot' .and. .not.pltxiBins(1)=='preview') then
    print *,"--Use given invalid option for pltxiBins(1)"
    flstopstatus = .true.
  endif
  if(.not.pltKLrealz(1)=='noplot' .and. .not.pltKLrealz(1)=='plot' .and. .not.pltKLrealz(1)=='preview') then
    print *,"--Use given invalid option for pltKLrealz(1)"
    flstopstatus = .true.
  endif
  if(.not.pltgenrealz(1)=='noplot' .and. .not.pltgenrealz(1)=='plot' .and. .not.pltgenrealz(1)=='preview') then
    print *,"--Use given invalid option for pltgenrealz(1)"
    flstopstatus = .true.
  endif
  if(.not.pltCo(1)=='noplot' .and. .not.pltCo(1)=='plot' .and. .not.pltCo(1)=='preview') then
    print *,"--Use given invalid option for pltCo(1)"
    flstopstatus = .true.
  endif
  if(.not.pltflux(1)=='noplot' .and. .not.pltflux(1)=='plot' .and. .not.pltflux(1)=='preview') then
    print *,"--Use given invalid option for pltflux(1)"
    flstopstatus = .true.
  endif
  if(.not.pltmatflux=='noplot' .and. .not.pltmatflux=='plot' .and. .not.pltmatflux=='preview') then
    print *,"--Use given invalid option for pltmatflux"
    flstopstatus = .true.
  endif


  if(flstopstatus) stop 'killed'
  if(flsleep) call sleep(4)

  call system("./saveinput.sh") !if pass tests, log input file for reference based on date/time
  call system("echo 'Input tests passed and run started:' $(date); echo")

  end subroutine read_test_inputstoc




  subroutine global_allocate
  !This subroutine allocates and initializes all global variables
  use rngvars, only: rngappnum, rngseed, mts
  use mt_stream, only: set_mt19937, new
  use genRealzvars, only: lam, P, slen, numRealz, numPath, sumPath, sqrPath, largesti, &
                          totLength, lamcs1, lamca1, lamcs2, lamca2, sig, aves1, &
                          avea1, aves2, avea2, scatrat, numPosRealz, numNegRealz, atmixscatrat, &
                          vars1, vara1, vars2, vara2, atmixsig, chgeomtype, &
                          GBlamcs1,GBlamca1,GBlamcs2,GBlamca2, &
                          GBaves1, GBavea1, GBaves2, GBavea2, GBvars1, GBvara1, GBvars2, GBvara2
  use KLvars, only: KLrnumpoints, pltKLrealznumof, chGausstype, Corropts, &
                    KLrxisig, alphas1, alphaa1, alphas2, alphaa2, Aks1, Aka1, Aks2, Aka2, &
                    Eigs1, Eiga1, Eigs2, Eiga2, lamctypes1, lamctypea1, lamctypes2, lamctypea2, &
                    pltCo, numEigss1, numEigsa1, numEigss2, numEigsa2, fls1, fla1, fls2, fla2, &
                    xi, xis1, xia1, xis2, xia2, corrinds1, corrinda1, &!corrinds2, corrinda2, &
                    pltKLrealzarray, pltKLrealz, numNystroms1, numNystroma1, numNystroms2, numNystroma2, &
                    eigvecss1, eigvecsa1, eigvecss2, eigvecsa2
  use MCvars, only: fluxfaces, numParts, stocMC_reflection, stocMC_transmission, &
                    stocMC_absorption, LPamnumParts, stocMC_fluxall, chTrantype, &
                    stocMC_fluxmat1, stocMC_fluxmat2, pltflux, pltmatflux, areapnsamp, &
                    fluxnumcells, flfluxplot, LPamMCsums, transmit, reflect, absorb, &
                    numpnSamp, radtrans_int, Wood_rej, flfluxplotall, flfluxplotmat, &
                    fluxall, numPartsperj, stocMC_reflectionPCE, stocMC_transmissionPCE, &
                    stocMC_fluxallPCE
  use UQvars, only: UQwgts, Qs, chUQtype, numPCEcoefs, PCEcoefsrefl, PCEcoefstran, PCEcoefscells, &
                    PCEcells, numPCEcells, PCEorder, flPCErefl, flPCEtran, numUQdims, samplePCExis, &
                    numPCEQoIsamps, numPCElocations, PCEreflsamples, PCEtransamples, PCEcellssamples
  use utilities, only: exponentialfit, nCr
  integer :: i

  !initialize rngvars
  rngappnum  = 0
  !initialize Mersenne Twister rngvars
  call set_mt19937
  call new(mts)

  !allocate and initialize genRealzvars
  if(chUQtype=='SC' .or. chUQtype=='PCE') numRealz = product(Qs)
  numPosRealz= 0
  numNegRealz= 0

  !allocate  KLresearch variables
  if(.not. fls1) numEigss1 = 0
  if(.not. fla1) numEigsa1 = 0
  if(.not. fls2) numEigss2 = 0
  if(.not. fla2) numEigsa2 = 0
  if(chTrantype=='KLWood' .or. chTrantype=='GaussKL' .or. &
     Corropts(1).ne.'noplot' .or. pltCo(1).ne.'noplot' .or. pltKLrealz(1).ne.'noplot') then
    if(fls1) allocate(Eigs1(numEigss1))
    if(fla1) allocate(Eiga1(numEigsa1))
    if(fls2) allocate(Eigs2(numEigss2))
    if(fla2) allocate(Eiga2(numEigsa2))
    if(fls1) allocate(eigvecss1(numEigss1,numNystroms1))
    if(fla1) allocate(eigvecsa1(numEigsa1,numNystroma1))
    if(fls2) allocate(eigvecss2(numEigss2,numNystroms2))
    if(fla2) allocate(eigvecsa2(numEigsa2,numNystroma2))
  endif

  if(chgeomtype=='contin') then  !Gauss-based input
    if(fls1) lamcs1       = GBlamcs1
    if(fla1) lamca1       = GBlamca1
    if(fls2) lamcs2       = GBlamcs2
    if(fla2) lamca2       = GBlamca2
    if(chGausstype=='Gaus') then
      if(fls1) aves1 = GBaves1
      if(fls1) vars1 = GBvars1

      if(fla1) avea1 = GBavea1
      if(fla1) vara1 = GBvara1

      if(fls2) aves2 = GBaves2
      if(fls2) vars2 = GBvars2

      if(fla2) avea2 = GBavea2
      if(fla2) vara2 = GBvara2
    elseif(chGausstype=='LogN') then
      if(fls1) aves1 = log( GBaves1**2 / sqrt( GBvars1 + GBaves1**2 ) )
      if(fls1) vars1 = log( GBvars1    / GBaves1**2 + 1d0 )

      if(fla1) avea1 = log( GBavea1**2 / sqrt( GBvara1 + GBavea1**2 ) )
      if(fla1) vara1 = log( GBvara1    / GBavea1**2 + 1d0 )

      if(fls2) aves2 = log( GBaves2**2 / sqrt( GBvars2 + GBaves2**2 ) )
      if(fls2) vars2 = log( GBvars2    / GBaves2**2 + 1d0 )

      if(fla2) avea2 = log( GBavea2**2 / sqrt( GBvara2 + GBavea2**2 ) )
      if(fla2) vara2 = log( GBvara2    / GBavea2**2 + 1d0 )

      if(fls1 .and. lamctypes1=='fitlamc') lamcs1 = exponentialfit(slen,1d0+GBvars1/GBaves1**2,lamcs1)
      if(fla1 .and. lamctypea1=='fitlamc') lamca1 = exponentialfit(slen,1d0+GBvara1/GBavea1**2,lamca1)
      if(fls2 .and. lamctypes2=='fitlamc') lamcs2 = exponentialfit(slen,1d0+GBvars2/GBaves2**2,lamcs2)
      if(fla2 .and. lamctypea2=='fitlamc') lamca2 = exponentialfit(slen,1d0+GBvara2/GBavea2**2,lamca2)
    endif

  elseif(chgeomtype=='binary') then
    numPath    = 0  !setup Markov material tallies
    sumPath    = 0d0
    sqrPath    = 0d0
    largesti   = 0
    totLength  = 0d0
    P(1)       = lam(1)/(lam(1)+lam(2)) !calc probabilities
    P(2)       = lam(2)/(lam(1)+lam(2))
    lamcs1     = (lam(1)*lam(2))/(lam(1)+lam(2))
    aves1 = P(1)*     scatrat(1) *sig(1) + P(2)*     scatrat(2) *sig(2)
    avea1 = P(1)*(1d0-scatrat(1))*sig(1) + P(2)*(1d0-scatrat(2))*sig(2)
    vars1 = P(1)*P(2) * (sig(1)*     scatrat(1)  - sig(2)*     scatrat(2)) **2
    vara1 = P(1)*P(2) * (sig(1)*(1d0-scatrat(1)) - sig(2)*(1d0-scatrat(2)))**2
    if(chTrantype=='KLWood') then
      if( (sig(1)*scatrat(1)-sig(2)*scatrat(2)>0d0 .and. sig(1)*(1d0-scatrat(1))-sig(2)*(1d0-scatrat(2))>0d0) .or. &
          (sig(1)*scatrat(1)-sig(2)*scatrat(2)<0d0 .and. sig(1)*(1d0-scatrat(1))-sig(2)*(1d0-scatrat(2))<0d0) ) then
        corrinds1 = 1
        corrinda1 = 1
      else
        corrinds1 = 1
        corrinda1 = -1
      endif
    endif
    if(chTrantype=='atmixMC') then
      atmixsig     =   P(1)*sig(1)            + P(2)*sig(2)
      atmixscatrat = ( P(1)*sig(1)*scatrat(1) + P(2)*sig(2)*scatrat(2) ) / atmixsig
    endif
  endif


  !allocate UQ variables
  allocate(UQwgts(numRealz))
  UQwgts = 0d0
  if(chUQtype=='PCE') then
    numPCEcoefs = nCr(PCEorder+numUQdims,PCEorder)
    if(flPCErefl) then
      allocate(PCEcoefsrefl(numPCEcoefs))
      PCEcoefsrefl = 0d0
    endif
    if(flPCEtran) then
      allocate(PCEcoefstran(numPCEcoefs))
      PCEcoefstran = 0d0
    endif
    numPCElocations = numPCEcells
    if(flPCErefl) numPCElocations = numPCElocations + 1
    if(flPCEtran) numPCElocations = numPCelocations + 1
    if(numPCEcells>0) then
      allocate(PCEcoefscells(numPCEcells,numPCEcoefs))
      PCEcoefscells = 0d0
      allocate(samplePCExis(numUQdims,numPCEQoIsamps,numPCElocations))
      samplePCExis = 0d0
    endif
    if(flPCErefl) then
      allocate(PCEreflsamples(numPCEQoIsamps))
      PCEreflsamples = 0d0
    endif
    if(flPCEtran) then
      allocate(PCEtransamples(numPCEQoIsamps))
      PCEtransamples = 0d0
    endif
    if(numPCEcells>0) then
      allocate(PCEcellssamples(numPCEcells,numPCEQoIsamps))
      PCEcellssamples = 0d0
    endif
  endif

  !allocate  KLresearch variables
  if(chTrantype=='KLWood' .or. chTrantype=='GaussKL' .or. &
     Corropts(1).ne.'noplot' .or. pltCo(1).ne.'noplot' .or. pltKLrealz(1).ne.'noplot') then
    allocate(alphas1(numEigss1))
    allocate(alphaa1(numEigsa1))
    allocate(alphas2(numEigss2))
    allocate(alphaa2(numEigsa2))
    allocate(Aks1(numEigss1))
    allocate(Aka1(numEigsa1))
    allocate(Aks2(numEigss2))
    allocate(Aka2(numEigsa2))
    allocate(xi(numRealz,max(numEigss1,numEigsa1,numEigss2,numEigsa2)))
  endif


  !allocate and initialize KLconstruction variables
  if(chTrantype=='KLWood' .or. chTrantype=='GaussKL' .or. pltKLrealz(1).ne.'noplot') then
    if(fls1) allocate(xis1(numRealz,numEigss1))
    if(fla1) allocate(xia1(numRealz,numEigsa1))
    if(fls2) allocate(xis2(numRealz,numEigss2))
    if(fla2) allocate(xia2(numRealz,numEigsa2))
    allocate(KLrxisig(KLrnumpoints))
    allocate(pltKLrealzarray(KLrnumpoints,pltKLrealznumof+1))
  endif



  !allocate/initialize MCvars
  if(.not.chTrantype=='None') then
    if(chTrantype=='LPMC' .or. chTrantype=='atmixMC') then
      numRealz = 1
      numParts = LPamnumParts
      if(.not.allocated(LPamMCsums)) allocate(LPamMCsums(3))
      LPamMCsums =0.0d0
    endif
    if(chTrantype=='radMC' .or. chTrantype=='radWood' .or. chTrantype=='KLWood' .or. chTrantype=='GaussKL') then
      allocate(numPartsperj(numRealz))
      numPartsperj = 0
    endif
    allocate(stocMC_reflection(3))   !global MC variables for each method
    allocate(stocMC_transmission(3)) !rank 3 holds 1=average, 2=deviation, 3=SEM
    allocate(stocMC_absorption(3))
    stocMC_reflection   = 0.0d0
    stocMC_transmission = 0.0d0
    stocMC_absorption   = 0.0d0
    if(chUQtype=='PCE') then
      allocate(stocMC_reflectionPCE(3))
      allocate(stocMC_transmissionPCE(3))
      stocMC_reflectionPCE   = 0d0
      stocMC_transmissionPCE = 0d0
    endif      
    if(.not.allocated(transmit)) allocate(transmit(numRealz)) !leakage tallies
    if(.not.allocated(reflect))  allocate(reflect(numRealz))
    if(.not.allocated(absorb))   allocate(absorb(numRealz))
    transmit     = 0.0d0
    reflect      = 0.0d0
    absorb       = 0.0d0

    if(.not.(pltflux(1)=='noplot' .and. pltmatflux=='noplot')) flfluxplot = .true.
    if(flfluxplot) then   !alloc and init flux cells
      allocate(fluxfaces(fluxnumcells+1))
      fluxfaces = 0.0d0
      do i=1,fluxnumcells+1
        fluxfaces(i) = (slen/fluxnumcells) * (i-1)
      enddo
    endif
    if(.not.pltflux(1)=='noplot' .or. &!mat irrespective flux allocations
      (flfluxplot .and. (chTrantype=='LPMC' .or. chTrantype=='atmixMC')) ) then
      allocate(stocMC_fluxall(fluxnumcells,2))
      stocMC_fluxall = 0.0d0
      if(chUQtype=='PCE' .and. numPCEcells>0) then
        allocate(stocMC_fluxallPCE(numPCEcells,3))
        stocMC_fluxallPCE = 0d0
      endif
    endif
    if(.not.pltmatflux=='noplot') then !mat respective flux allocations
      allocate(stocMC_fluxmat1(fluxnumcells,2))
      allocate(stocMC_fluxmat2(fluxnumcells,2))
      stocMC_fluxmat1 = 0.0d0
      stocMC_fluxmat2 = 0.0d0
    endif
    if(chTrantype=='radWood' .or. chTrantype=='KLWood' .or. & !rejection tally allocations
       chTrantype=='GaussKL'                                ) Wood_rej = 0
    radtrans_int = 0
    if(chTrantype=='KLWood' .or.  chTrantype=='GaussKL' ) then !negative xs transport tally allocations
      numPosRealz=0
      numNegRealz=0
      numpnSamp  =0
      areapnSamp =0.0d0
    endif
    if(.not. pltflux(1)=='noplot') flfluxplotall = .true.
    if(flfluxplotall) then
      allocate(fluxall(fluxnumcells,numRealz))
      fluxall = 0d0
    endif
    if(.not. pltmatflux=='noplot') flfluxplotmat = .true.

  endif

  end subroutine global_allocate




  subroutine clearreports
  call system("rm texts/*.out")
  end subroutine clearreports



  subroutine finalreport
  !compile and display on screen a final report of some data
  call system("cat texts/geominput.out > texts/finalreport.out")
  call system("test -e texts/Woodnegstats.out && cat texts/Woodnegstats.out >> texts/finalreport.out")
  call system("test -e texts/MCleakage.out && cat texts/MCleakage.out >> texts/finalreport.out")
  call system("cat texts/finalreport.out")
  end subroutine finalreport



  subroutine GBcase_load
  !load special cases; not take input from input file.
  !right now the only special cases are Fichtl 1 and Fichtl 2
  use genRealzvars, only: GBaves1, GBavea1, GBaves2, GBavea2, GBvars1, GBvara1, GBvars2, GBvara2, &
                          GBlamcs1, GBlamca1, GBlamcs2, GBlamca2, slen
  use KLvars, only: chGBcase, corrinds1, corrinda1, corrinds2, corrinda2, numEigss1, numEigsa1, numEigss2, numEigsa2, &
                    fls1, fla1, fls2, fla2
  use MCvars, only: sourceType

  if(chGBcase=='ANS16l' .or. chGBcase=='ANS16s') then
    fls1 = .true.
    fla1 = .true.
    fls2 = .false.
    fla2 = .false.

    GBavea1   = 0.75d0  
    GBaves1   = 0.75d0
    if(chGBcase=='ANS16l') then
      GBvara1   = 2.25d0
      GBvars1   = 2.25d0
    elseif(chGBcase=='ANS16s') then
      GBvara1   = 0.075d0
      GBvars1   = 0.075d0
    endif

    GBlamcs1  = 1.5d0
    GBlamca1  = 1.5d0
    slen      = 5.0d0

    numEigss1 = 4
    numEigsa1 = 4
    numEigss2 = 0
    numEigsa2 = 0

    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 0
    corrinda2 = 0

    sourceType  = 'leftbeam'
  endif
  

  if(chGBcase=='f1' .or. chGBcase=='f2') then
    fls1 = .true.
    fla1 = .true.
    fls2 = .false.
    fla2 = .false.

    if(chGBcase=='f1') then
      GBavea1    = 2.5d0  
      GBvara1    = 0.5d0
      GBaves1    = 2.5d0
      GBvars1    = 0.5d0
    elseif(chGBcase=='f2') then
      GBavea1    = 0.5d0  
      GBvara1    = 0.02d0
      GBaves1    = 4.5d0
      GBvars1    = 1.62d0
    endif

    GBlamcs1  = 1.0d0
    GBlamca1  = 1.0d0
    slen      = 5.0d0

    numEigss1 = 5
    numEigsa1 = 5
    numEigss2 = 0
    numEigsa2 = 0

    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0

    sourceType  = 'leftbeam'
  endif

  if(chGBcase=='A1.1') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0.90909090909090995d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892562704d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 0.10000000000000001d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A1.1ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0.90909090909090995d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892562704d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 0.10000000000000001d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A1.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90909090909090995d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 7.43801652892562704d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A1.2') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0.90909090909090995d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892562704d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 1.00000000000000000d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A1.2ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0.90909090909090995d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892562704d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 1.00000000000000000d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A1.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90909090909090995d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 7.43801652892562704d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A1.3') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0.90909090909090995d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892562704d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 10.00000000000000000d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A1.3ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0.90909090909090995d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892562704d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 10.00000000000000000d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A1.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90909090909090995d0
    GBavea1   = 0.09090909090909090d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 7.43801652892562704d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A2.1') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea2   = 0.90909090909090995d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892562704d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A2.1ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea2   = 0.90909090909090995d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892562704d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A2.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea1   = 0.90909090909090995d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara1   = 7.43801652892562704d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A2.2') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea2   = 0.90909090909090995d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892562704d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A2.2ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea2   = 0.90909090909090995d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892562704d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A2.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea1   = 0.90909090909090995d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara1   = 7.43801652892562704d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A2.3') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea2   = 0.90909090909090995d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892562704d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A2.3ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea2   = 0.90909090909090995d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892562704d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A2.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.09090909090909090d0
    GBavea1   = 0.90909090909090995d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara1   = 7.43801652892562704d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A3.1') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818180d0
    GBavea2   = 0.09090909090909098d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181901d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925623d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975931d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A3.1ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818180d0
    GBavea2   = 0.09090909090909098d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181901d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925623d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975931d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A3.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90000000000000080d0
    GBavea1   = 0.10000000000000006d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 5.89165289256198932d0
    GBvara1   = 0.07273645546372823d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A3.2') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818180d0
    GBavea2   = 0.09090909090909098d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181901d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925623d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975931d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A3.2ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818180d0
    GBavea2   = 0.09090909090909098d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181901d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925623d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975931d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A3.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90000000000000080d0
    GBavea1   = 0.10000000000000006d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 5.89165289256198932d0
    GBvara1   = 0.07273645546372823d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A3.3') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818180d0
    GBavea2   = 0.09090909090909098d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181901d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925623d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975931d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A3.3ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818180d0
    GBavea2   = 0.09090909090909098d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181901d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925623d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975931d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca2  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0.09899999999999999d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A3.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90000000000000080d0
    GBavea1   = 0.10000000000000006d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 5.89165289256198932d0
    GBvara1   = 0.07273645546372823d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.09899999999999999d0
    GBlamca1  = 0.09899999999999999d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A4.1') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0.90909090909090895d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892561993d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 0.10000000000000001d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A4.1ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0.90909090909090895d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892561993d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 0.10000000000000001d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A4.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90909090909090895d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 7.43801652892561993d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A4.2') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0.90909090909090895d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892561993d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 1.00000000000000000d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A4.2ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0.90909090909090895d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892561993d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 1.00000000000000000d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A4.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90909090909090895d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 7.43801652892561993d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A4.3') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0.90909090909090895d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892561993d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 10.00000000000000000d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A4.3ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0.90909090909090895d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 7.43801652892561993d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 10.00000000000000000d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A4.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90909090909090895d0
    GBavea1   = 0.09090909090909091d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 7.43801652892561993d0
    GBvara1   = 0.00091827364554637d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A5.1') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea2   = 0.90909090909090895d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892561993d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A5.1ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea2   = 0.90909090909090895d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892561993d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A5.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea1   = 0.90909090909090895d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara1   = 7.43801652892561993d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A5.2') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea2   = 0.90909090909090895d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892561993d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A5.2ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea2   = 0.90909090909090895d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892561993d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A5.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea1   = 0.90909090909090895d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara1   = 7.43801652892561993d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A5.3') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea2   = 0.90909090909090895d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892561993d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A5.3ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea2   = 0.90909090909090895d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara2   = 7.43801652892561993d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A5.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.09090909090909091d0
    GBavea1   = 0.90909090909090895d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00091827364554637d0
    GBvara1   = 7.43801652892561993d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A6.1') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818182d0
    GBavea2   = 0.09090909090909087d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181812d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925616d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975398d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A6.1ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818182d0
    GBavea2   = 0.09090909090909087d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181812d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925616d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975398d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A6.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.89999999999999991d0
    GBavea1   = 0.09999999999999996d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 5.89165289256198399d0
    GBvara1   = 0.07273645546372816d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A6.2') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818182d0
    GBavea2   = 0.09090909090909087d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181812d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925616d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975398d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A6.2ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818182d0
    GBavea2   = 0.09090909090909087d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181812d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925616d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975398d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A6.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.89999999999999991d0
    GBavea1   = 0.09999999999999996d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 5.89165289256198399d0
    GBvara1   = 0.07273645546372816d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A6.3') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818182d0
    GBavea2   = 0.09090909090909087d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181812d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925616d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975398d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A6.3ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.08181818181818182d0
    GBavea2   = 0.09090909090909087d0
    GBavea1   = 0.00909090909090909d0
    GBaves2   = 0.81818181818181812d0
    GBvars1   = 0.00074380165289256d0
    GBvara2   = 0.07438016528925616d0
    GBvara1   = 0.00000918273645546d0
    GBvars2   = 6.02479338842975398d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca2  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0.99000000000000010d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A6.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.89999999999999991d0
    GBavea1   = 0.09999999999999996d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 5.89165289256198399d0
    GBvara1   = 0.07273645546372816d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 0.99000000000000010d0
    GBlamca1  = 0.99000000000000010d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A7.1') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0.99009900990099009d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0.98029604940692083d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 0.10000000000000001d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A7.1ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0.99009900990099009d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0.98029604940692083d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 0.10000000000000001d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A7.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.99009900990099009d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.98029604940692083d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A7.2') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0.99009900990099009d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0.98029604940692083d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 1.00000000000000000d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A7.2ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0.99009900990099009d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0.98029604940692083d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 1.00000000000000000d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A7.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.99009900990099009d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.98029604940692083d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A7.3') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0.99009900990099009d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0.98029604940692083d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 10.00000000000000000d0
    corrinds1 = 0
    corrinda2 = 0
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A7.3ind') then
    sourceType= 'leftiso'
    fls1      = .false.
    fla2      = .false.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0d0
    GBavea2   = 0d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0.99009900990099009d0
    GBvars1   = 0d0
    GBvara2   = 0d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0.98029604940692083d0
    GBlamcs1  = 0d0
    GBlamca2  = 0d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 10.00000000000000000d0
    corrinda1 = 1
    corrinds2 = 2
    corrinds1 = 0
    corrinda2 = 0
  elseif(chGBcase=='A7.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.99009900990099009d0
    GBavea1   = 0.00990099009900990d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.98029604940692083d0
    GBvara1   = 0.00009802960494069d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A8.1') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea2   = 0.99009900990099009d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara2   = 0.98029604940692083d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A8.1ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea2   = 0.99009900990099009d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara2   = 0.98029604940692083d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A8.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea1   = 0.99009900990099009d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara1   = 0.98029604940692083d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A8.2') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea2   = 0.99009900990099009d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara2   = 0.98029604940692083d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A8.2ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea2   = 0.99009900990099009d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara2   = 0.98029604940692083d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A8.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea1   = 0.99009900990099009d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara1   = 0.98029604940692083d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A8.3') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea2   = 0.99009900990099009d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara2   = 0.98029604940692083d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A8.3ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .false.
    fls2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea2   = 0.99009900990099009d0
    GBavea1   = 0d0
    GBaves2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara2   = 0.98029604940692083d0
    GBvara1   = 0d0
    GBvars2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 0d0
    GBlamcs2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = 2
    corrinda1 = 0
    corrinds2 = 0
  elseif(chGBcase=='A8.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.00990099009900990d0
    GBavea1   = 0.99009900990099009d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.00009802960494069d0
    GBvara1   = 0.98029604940692083d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = -1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A9.1') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.00891089108910891d0
    GBavea2   = 0.09900990099009899d0
    GBavea1   = 0.00099009900990099d0
    GBaves2   = 0.89108910891089110d0
    GBvars1   = 0.00007940398000196d0
    GBvara2   = 0.00980296049406921d0
    GBvara1   = 0.00000098029604941d0
    GBvars2   = 0.79403980001960595d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A9.1ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.00891089108910891d0
    GBavea2   = 0.09900990099009899d0
    GBavea1   = 0.00099009900990099d0
    GBaves2   = 0.89108910891089110d0
    GBvars1   = 0.00007940398000196d0
    GBvara2   = 0.00980296049406921d0
    GBvara1   = 0.00000098029604941d0
    GBvars2   = 0.79403980001960595d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A9.1two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90000000000000002d0
    GBavea1   = 0.09999999999999998d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.77823840799921573d0
    GBvara1   = 0.00960788158023723d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 0.10000000000000001d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A9.2') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.00891089108910891d0
    GBavea2   = 0.09900990099009899d0
    GBavea1   = 0.00099009900990099d0
    GBaves2   = 0.89108910891089110d0
    GBvars1   = 0.00007940398000196d0
    GBvara2   = 0.00980296049406921d0
    GBvara1   = 0.00000098029604941d0
    GBvars2   = 0.79403980001960595d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A9.2ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.00891089108910891d0
    GBavea2   = 0.09900990099009899d0
    GBavea1   = 0.00099009900990099d0
    GBaves2   = 0.89108910891089110d0
    GBvars1   = 0.00007940398000196d0
    GBvara2   = 0.00980296049406921d0
    GBvara1   = 0.00000098029604941d0
    GBvars2   = 0.79403980001960595d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A9.2two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90000000000000002d0
    GBavea1   = 0.09999999999999998d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.77823840799921573d0
    GBvara1   = 0.00960788158023723d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 1.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  elseif(chGBcase=='A9.3') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.00891089108910891d0
    GBavea2   = 0.09900990099009899d0
    GBavea1   = 0.00099009900990099d0
    GBaves2   = 0.89108910891089110d0
    GBvars1   = 0.00007940398000196d0
    GBvara2   = 0.00980296049406921d0
    GBvara1   = 0.00000098029604941d0
    GBvars2   = 0.79403980001960595d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda2 = -1
    corrinda1 = 1
    corrinds2 = -1
  elseif(chGBcase=='A9.3ind') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla2      = .true.
    fla1      = .true.
    fls2      = .true.
    GBaves1   = 0.00891089108910891d0
    GBavea2   = 0.09900990099009899d0
    GBavea1   = 0.00099009900990099d0
    GBaves2   = 0.89108910891089110d0
    GBvars1   = 0.00007940398000196d0
    GBvara2   = 0.00980296049406921d0
    GBvara1   = 0.00000098029604941d0
    GBvars2   = 0.79403980001960595d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca2  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 2.52499999999999991d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 2
    corrinds2 = 3
    corrinda2 = 4
  elseif(chGBcase=='A9.3two') then
    sourceType= 'leftiso'
    fls1      = .true.
    fla1      = .true.
    fls2      = .false.
    fla2      = .false.
    GBaves1   = 0.90000000000000002d0
    GBavea1   = 0.09999999999999998d0
    GBaves2   = 0d0
    GBavea2   = 0d0
    GBvars1   = 0.77823840799921573d0
    GBvara1   = 0.00960788158023723d0
    GBvars2   = 0d0
    GBvara2   = 0d0
    GBlamcs1  = 2.52499999999999991d0
    GBlamca1  = 2.52499999999999991d0
    GBlamcs2  = 0d0
    GBlamca2  = 0d0
    slen      = 10.00000000000000000d0
    corrinds1 = 1
    corrinda1 = 1
    corrinds2 = 0
    corrinda2 = 0
  endif
print *,"chGBcase:",chGBcase
print *,fls1,fla1,fls2,fla2,GBaves1,GBavea1,GBaves2,GBavea2,GBvars1,GBvara1,GBvars2,GBvara2
print *,GBlamcs1,GBlamca1,GBlamcs2,GBlamca2,slen,corrinds1,corrinda1,corrinds2,corrinda2

  end subroutine GBcase_load



  subroutine Acase_load
  use genRealzvars, only:   Adamscase, sig, lam, scatrat, slen  
  use MCvars, only: rodOrplanar, ABreflection, ABtransmission, sourceType
  real(8) :: eps = 0.0001d0 !local var
    ! Loading specials cases from Adams and Pomraning Paper, 27 total cases
    sourceType = 'leftiso'
    if( Adamscase /= 0 ) then
      ! load sig and lam
      if( int(Adamscase)==1 .OR. int(Adamscase)==2 .OR. int(Adamscase)==3 ) then
        sig(1)=0.10101010d0
        sig(2)=9.0909091d0
        lam(1)=0.99d0
        lam(2)=0.11d0
      elseif( int(Adamscase)==4 .OR. int(Adamscase)==5 .OR. int(Adamscase)==6 ) then
        sig(1)=0.10101010d0
        sig(2)=9.0909091d0
        lam(1)=9.9d0
        lam(2)=1.1d0
      elseif( int(Adamscase)==7 .OR. int(Adamscase)==8 .OR. int(Adamscase)==9 ) then
        sig(1)=0.019801980d0
        sig(2)=1.9801980d0
        lam(1)=5.05d0
        lam(2)=5.05d0
      endif
      ! load scatrat
      if( int(Adamscase)==1 .OR. int(Adamscase)==4 .OR. int(Adamscase)==7 ) then
        scatrat(1)=0
        scatrat(2)=1
      elseif( int(Adamscase)==2 .OR. int(Adamscase)==5 .OR. int(Adamscase)==8 ) then
        scatrat(1)=1
        scatrat(2)=0
      elseif( int(Adamscase)==3 .OR. int(Adamscase)==6 .OR. int(Adamscase)==9 ) then
        scatrat(1)=0.9d0
        scatrat(2)=0.9d0
      endif
      ! load s
      if( int((Adamscase*10+0.005d0-(int(Adamscase))*10))==1 ) then
        slen=0.1d0
      elseif( int((Adamscase*10+0.005d0-(int(Adamscase))*10))==2 ) then
        slen=1.0d0
      elseif( int((Adamscase*10+0.005d0-(int(Adamscase))*10))==3 ) then
        slen=10.0d0
      endif
    endif

    ! Loading transmission and reflection results from Adams and Brantley Papers
    allocate(ABreflection(2,5))
    allocate(ABtransmission(2,5))
    ABreflection   = 0d0
    ABtransmission = 0d0
    if(rodOrplanar=='rod') then
      if( Adamscase<1.1d0+eps .and. Adamscase>1.1d0-eps ) then
        ABreflection(1,1)   = 0.0344d0 !AdMCave
        ABreflection(2,1)   = 0.0843d0 !AdMCdev
        ABreflection(1,2)   = 0.0337d0 !AdLPave

        ABtransmission(1,1) = 0.9566d0 !AdMCave
        ABtransmission(2,1) = 0.0818d0 !AdMCdev
        ABtransmission(1,2) = 0.9572d0 !AdLPave
      elseif( Adamscase<1.2d0+eps .and. Adamscase>1.2d0-eps ) then
        ABreflection(1,1)   = 0.2130d0 !AdMCave
        ABreflection(2,1)   = 0.2100d0 !AdMCdev
        ABreflection(1,2)   = 0.1948d0 !AdLPave

        ABtransmission(1,1) = 0.7006d0 !AdMCave
        ABtransmission(2,1) = 0.1982d0 !AdMCdev
        ABtransmission(1,2) = 0.7187d0 !AdLPave
      elseif( Adamscase<1.3d0+eps .and. Adamscase>1.3d0-eps ) then
        ABreflection(1,1)   = 0.4886d0 !AdMCave
        ABreflection(2,1)   = 0.1313d0 !AdMCdev
        ABreflection(1,2)   = 0.4393d0 !AdLPave

        ABtransmission(1,1) = 0.0526d0 !AdMCave
        ABtransmission(2,1) = 0.0339d0 !AdMCdev
        ABtransmission(1,2) = 0.0772d0 !AdLPave
      elseif( Adamscase<2.1d0+eps .and. Adamscase>2.1d0-eps ) then
        ABreflection(1,1)   = 0.0044d0 !AdMCave
        ABreflection(2,1)   = 0.0014d0 !AdMCdev
        ABreflection(1,2)   = 0.0044d0 !AdLPave

        ABtransmission(1,1) = 0.9285d0 !AdMCave
        ABtransmission(2,1) = 0.1618d0 !AdMCdev
        ABtransmission(1,2) = 0.9293d0 !AdLPave
      elseif( Adamscase<2.2d0+eps .and. Adamscase>2.2d0-eps ) then
        ABreflection(1,1)   = 0.0325d0 !AdMCave
        ABreflection(2,1)   = 0.01959d0 !AdMCdev
        ABreflection(1,2)   = 0.0292d0 !AdLPave

        ABtransmission(1,1) = 0.5699d0 !AdMCave
        ABtransmission(2,1) = 0.3494d0 !AdMCdev
        ABtransmission(1,2) = 0.5699d0 !AdLPave
      elseif( Adamscase<2.3d0+eps .and. Adamscase>2.3d0-eps ) then
        ABreflection(1,1)   = 0.0615d0 !AdMCave
        ABreflection(2,1)   = 0.0486d0 !AdMCdev
        ABreflection(1,2)   = 0.0443d0 !AdLPave

        ABtransmission(1,1) = 0.0047d0 !AdMCave
        ABtransmission(2,1) = 0.0223d0 !AdMCdev
        ABtransmission(1,2) = 0.0045d0 !AdLPave
      elseif( Adamscase<3.1d0+eps .and. Adamscase>3.1d0-eps ) then
        ABreflection(1,1)   = 0.0339d0 !AdMCave
        ABreflection(2,1)   = 0.0717d0 !AdMCdev
        ABreflection(1,2)   = 0.0333d0 !AdLPave

        ABtransmission(1,1) = 0.9563d0 !AdMCave
        ABtransmission(2,1) = 0.0938d0 !AdMCdev
        ABtransmission(1,2) = 0.9570d0 !AdLPave
      elseif( Adamscase<3.2d0+eps .and. Adamscase>3.2d0-eps ) then
        ABreflection(1,1)   = 0.2120d0 !AdMCave
        ABreflection(2,1)   = 0.1586d0 !AdMCdev
        ABreflection(1,2)   = 0.1929d0 !AdLPave

        ABtransmission(1,1) = 0.7017d0 !AdMCave
        ABtransmission(2,1) = 0.2437d0 !AdMCdev
        ABtransmission(1,2) = 0.7182d0 !AdLPave
      elseif( Adamscase<3.3d0+eps .and. Adamscase>3.3d0-eps ) then
        ABreflection(1,1)   = 0.5146d0 !AdMCave
        ABreflection(2,1)   = 0.0123d0 !AdMCdev
        ABreflection(1,2)   = 0.4338d0 !AdLPave

        ABtransmission(1,1) = 0.0557d0 !AdMCave
        ABtransmission(2,1) = 0.0623d0 !AdMCdev
        ABtransmission(1,2) = 0.0759d0 !AdLPave
      elseif( Adamscase<4.1d0+eps .and. Adamscase>4.1d0-eps ) then
        ABreflection(1,1)   = 0.0319d0 !AdMCave
        ABreflection(2,1)   = 0.0931d0 !AdMCdev
        ABreflection(1,2)   = 0.0315d0 !AdLPave

        ABtransmission(1,1) = 0.9591d0 !AdMCave
        ABtransmission(2,1) = 0.0901d0 !AdMCdev
        ABtransmission(1,2) = 0.9594d0 !AdLPave
      elseif( Adamscase<4.2d0+eps .and. Adamscase>4.2d0-eps ) then
        ABreflection(1,1)   = 0.1130d0 !AdMCave
        ABreflection(2,1)   = 0.2531d0 !AdMCdev
        ABreflection(1,2)   = 0.1003d0 !AdLPave

        ABtransmission(1,1) = 0.8006d0 !AdMCave
        ABtransmission(2,1) = 0.2295d0 !AdMCdev
        ABtransmission(1,2) = 0.8132d0 !AdLPave
      elseif( Adamscase<4.3d0+eps .and. Adamscase>4.3d0-eps ) then
        ABreflection(1,1)   = 0.2841d0 !AdMCave
        ABreflection(2,1)   = 0.2876d0 !AdMCdev
        ABreflection(1,2)   = 0.2158d0 !AdLPave

        ABtransmission(1,1) = 0.1774d0 !AdMCave
        ABtransmission(2,1) = 0.1437d0 !AdMCdev
        ABtransmission(1,2) = 0.2294d0 !AdLPave
      elseif( Adamscase<5.1d0+eps .and. Adamscase>5.1d0-eps ) then
        ABreflection(1,1)   = 0.0045d0 !AdMCave
        ABreflection(2,1)   = 0.0015d0 !AdMCdev
        ABreflection(1,2)   = 0.0045d0 !AdLPave

        ABtransmission(1,1) = 0.9343d0 !AdMCave
        ABtransmission(2,1) = 0.1768d0 !AdMCdev
        ABtransmission(1,2) = 0.9350d0 !AdLPave
      elseif( Adamscase<5.2d0+eps .and. Adamscase>5.2d0-eps ) then
        ABreflection(1,1)   = 0.0414d0 !AdMCave
        ABreflection(2,1)   = 0.0154d0 !AdMCdev
        ABreflection(1,2)   = 0.0402d0 !AdLPave

        ABtransmission(1,1) = 0.7948d0 !AdMCave
        ABtransmission(2,1) = 0.3384d0 !AdMCdev
        ABtransmission(1,2) = 0.7970d0 !AdLPave
      elseif( Adamscase<5.3d0+eps .and. Adamscase>5.3d0-eps ) then
        ABreflection(1,1)   = 0.2078d0 !AdMCave
        ABreflection(2,1)   = 0.1243d0 !AdMCdev
        ABreflection(1,2)   = 0.1568d0 !AdLPave

        ABtransmission(1,1) = 0.2430d0 !AdMCave
        ABtransmission(2,1) = 0.3067d0 !AdMCdev
        ABtransmission(1,2) = 0.2352d0 !AdLPave
      elseif( Adamscase<6.1d0+eps .and. Adamscase>6.1d0-eps ) then
        ABreflection(1,1)   = 0.0315d0 !AdMCave
        ABreflection(2,1)   = 0.0785d0 !AdMCdev
        ABreflection(1,2)   = 0.0311d0 !AdLPave

        ABtransmission(1,1) = 0.9589d0 !AdMCave
        ABtransmission(2,1) = 0.1039d0 !AdMCdev
        ABtransmission(1,2) = 0.9593d0 !AdLPave
      elseif( Adamscase<6.2d0+eps .and. Adamscase>6.2d0-eps ) then
        ABreflection(1,1)   = 0.1181d0 !AdMCave
        ABreflection(2,1)   = 0.1632d0 !AdMCdev
        ABreflection(1,2)   = 0.1045d0 !AdLPave

        ABtransmission(1,1) = 0.8181d0 !AdMCave
        ABtransmission(2,1) = 0.2873d0 !AdMCdev
        ABtransmission(1,2) = 0.8270d0 !AdLPave
      elseif( Adamscase<6.3d0+eps .and. Adamscase>6.3d0-eps ) then
        ABreflection(1,1)   = 0.4301d0 !AdMCave
        ABreflection(2,1)   = 0.1066d0 !AdMCdev
        ABreflection(1,2)   = 0.2883d0 !AdLPave

        ABtransmission(1,1) = 0.2658d0 !AdMCave
        ABtransmission(2,1) = 0.2728d0 !AdMCdev
        ABtransmission(1,2) = 0.2961d0 !AdLPave
      elseif( Adamscase<7.1d0+eps .and. Adamscase>7.1d0-eps ) then
        ABreflection(1,1)   = 0.0450d0 !AdMCave
        ABreflection(2,1)   = 0.0448d0 !AdMCdev
        ABreflection(1,2)   = 0.0451d0 !AdLPave

        ABtransmission(1,1) = 0.9540d0 !AdMCave
        ABtransmission(2,1) = 0.0438d0 !AdMCdev
        ABtransmission(1,2) = 0.9539d0 !AdLPave
      elseif( Adamscase<7.2d0+eps .and. Adamscase>7.2d0-eps ) then
        ABreflection(1,1)   = 0.2586d0 !AdMCave
        ABreflection(2,1)   = 0.2338d0 !AdMCdev
        ABreflection(1,2)   = 0.2562d0 !AdLPave

        ABtransmission(1,1) = 0.7316d0 !AdMCave
        ABtransmission(2,1) = 0.2247d0 !AdMCdev
        ABtransmission(1,2) = 0.7340d0 !AdLPave
      elseif( Adamscase<7.3d0+eps .and. Adamscase>7.3d0-eps ) then
        ABreflection(1,1)   = 0.6804d0 !AdMCave
        ABreflection(2,1)   = 0.2528d0 !AdMCdev
        ABreflection(1,2)   = 0.6034d0 !AdLPave

        ABtransmission(1,1) = 0.2328d0 !AdMCave
        ABtransmission(2,1) = 0.2029d0 !AdMCdev
        ABtransmission(1,2) = 0.3063d0 !AdLPave
      elseif( Adamscase<8.1d0+eps .and. Adamscase>8.1d0-eps ) then
        ABreflection(1,1)   = 0.0005d0 !AdMCave
        ABreflection(2,1)   = 0.0005d0 !AdMCdev
        ABreflection(1,2)   = 0.0005d0 !AdLPave

        ABtransmission(1,1) = 0.9097d0 !AdMCave
        ABtransmission(2,1) = 0.0888d0 !AdMCdev
        ABtransmission(1,2) = 0.9093d0 !AdLPave
      elseif( Adamscase<8.2d0+eps .and. Adamscase>8.2d0-eps ) then
        ABreflection(1,1)   = 0.0046d0 !AdMCave
        ABreflection(2,1)   = 0.0047d0 !AdMCdev
        ABreflection(1,2)   = 0.0045d0 !AdLPave

        ABtransmission(1,1) = 0.5396d0 !AdMCave
        ABtransmission(2,1) = 0.4027d0 !AdMCdev
        ABtransmission(1,2) = 0.5396d0 !AdLPave
      elseif( Adamscase<8.3d0+eps .and. Adamscase>8.3d0-eps ) then
        ABreflection(1,1)   = 0.0217d0 !AdMCave
        ABreflection(2,1)   = 0.0294d0 !AdMCdev
        ABreflection(1,2)   = 0.0148d0 !AdLPave

        ABtransmission(1,1) = 0.0910d0 !AdMCave
        ABtransmission(2,1) = 0.2467d0 !AdMCdev
        ABtransmission(1,2) = 0.0913d0 !AdLPave
      elseif( Adamscase<9.1d0+eps .and. Adamscase>9.1d0-eps ) then
        ABreflection(1,1)   = 0.0406d0 !AdMCave
        ABreflection(2,1)   = 0.0394d0 !AdMCdev
        ABreflection(1,2)   = 0.0406d0 !AdLPave

        ABtransmission(1,1) = 0.9495d0 !AdMCave
        ABtransmission(2,1) = 0.0491d0 !AdMCdev
        ABtransmission(1,2) = 0.9495d0 !AdLPave
      elseif( Adamscase<9.2d0+eps .and. Adamscase>9.2d0-eps ) then
        ABreflection(1,1)   = 0.2152d0 !AdMCave
        ABreflection(2,1)   = 0.1852d0 !AdMCdev
        ABreflection(1,2)   = 0.2129d0 !AdLPave

        ABtransmission(1,1) = 0.6955d0 !AdMCave
        ABtransmission(2,1) = 0.2661d0 !AdMCdev
        ABtransmission(1,2) = 0.6976d0 !AdLPave
      elseif( Adamscase<9.3d0+eps .and. Adamscase>9.3d0-eps ) then
        ABreflection(1,1)   = 0.4688d0 !AdMCave
        ABreflection(2,1)   = 0.1214d0 !AdMCdev
        ABreflection(1,2)   = 0.3693d0 !AdLPave

        ABtransmission(1,1) = 0.1510d0 !AdMCave
        ABtransmission(2,1) = 0.2566d0 !AdMCdev
        ABtransmission(1,2) = 0.1798d0 !AdLPave
      endif
    elseif(rodOrplanar=='planar') then
      if( Adamscase<1.1d0+eps .and. Adamscase>1.1d0-eps ) then
        ABreflection(1,1)   = 0.0491d0 !AdMCave
        ABreflection(2,1)   = 0.1182d0 !AdMCdev
        ABreflection(1,2)   = 0.0479d0 !AdLPave
        ABreflection(1,3)   = 0.04874d0 !BrMCave
        ABreflection(1,4)   = 0.04768d0 !BrLPave
        ABreflection(1,5)   = 0.07544d0 !Bratmixave

        ABtransmission(1,1) = 0.9331d0 !AdMCave
        ABtransmission(2,1) = 0.1132d0 !AdMCdev
        ABtransmission(1,2) = 0.9343d0 !AdLPave
        ABtransmission(1,3) = 0.93369d0 !BrMCave
        ABtransmission(1,4) = 0.93463d0 !BrLPave
        ABtransmission(1,5) = 0.90690d0 !Bratmixave
      elseif( Adamscase<1.2d0+eps .and. Adamscase>1.2d0-eps ) then
        ABreflection(1,1)   = 0.2495d0 !AdMCave
        ABreflection(2,1)   = 0.2354d0 !AdMCdev
        ABreflection(1,2)   = 0.2187d0 !AdLPave
        ABreflection(1,3)   = 0.25069d0 !BrMCave
        ABreflection(1,4)   = 0.21888d0 !BrLPave
        ABreflection(1,5)   = 0.36040d0 !Bratmixave

        ABtransmission(1,1) = 0.5950d0 !AdMCave
        ABtransmission(2,1) = 0.2143d0 !AdMCdev
        ABtransmission(1,2) = 0.6254d0 !AdLPave
        ABtransmission(1,3) = 0.59607d0 !BrMCave
        ABtransmission(1,4) = 0.62741d0 !BrLPave
        ABtransmission(1,5) = 0.48096d0 !Bratmixave
      elseif( Adamscase<1.3d0+eps .and. Adamscase>1.3d0-eps ) then
        ABreflection(1,1)   = 0.4342d0 !AdMCave
        ABreflection(2,1)   = 0.1616d0 !AdMCdev
        ABreflection(1,2)   = 0.3760d0 !AdLPave
        ABreflection(1,3)   = 0.43634d0 !BrMCave
        ABreflection(1,4)   = 0.37815d0 !BrLPave
        ABreflection(1,5)   = 0.49564d0 !Bratmixave

        ABtransmission(1,1) = 0.0146d0 !AdMCave
        ABtransmission(2,1) = 0.0152d0 !AdMCdev
        ABtransmission(1,2) = 0.0259d0 !AdLPave
        ABtransmission(1,3) = 0.01486d0 !BrMCave
        ABtransmission(1,4) = 0.02642d0 !BrLPave
        ABtransmission(1,5) = 0.00474d0 !Bratmixave
      elseif( Adamscase<2.1d0+eps .and. Adamscase>2.1d0-eps ) then
        ABreflection(1,1)   = 0.0087d0 !AdMCave
        ABreflection(2,1)   = 0.0029d0 !AdMCdev
        ABreflection(1,2)   = 0.0086d0 !AdLPave
        ABreflection(1,3)   = 0.00862d0 !BrMCave
        ABreflection(1,4)   = 0.00849d0 !BrLPave
        ABreflection(1,5)   = 0.00654d0 !Bratmixave

        ABtransmission(1,1) = 0.9014d0 !AdMCave
        ABtransmission(2,1) = 0.2111d0 !AdMCdev
        ABtransmission(1,2) = 0.9004d0 !AdLPave
        ABtransmission(1,3) = 0.90118d0 !BrMCave
        ABtransmission(1,4) = 0.90080d0 !BrLPave
        ABtransmission(1,5) = 0.83900d0 !Bratmixave
      elseif( Adamscase<2.2d0+eps .and. Adamscase>2.2d0-eps ) then
        ABreflection(1,1)   = 0.0548d0 !AdMCave
        ABreflection(2,1)   = 0.0307d0 !AdMCdev
        ABreflection(1,2)   = 0.0460d0 !AdLPave
        ABreflection(1,3)   = 0.05418d0 !BrMCave
        ABreflection(1,4)   = 0.04552d0 !BrLPave
        ABreflection(1,5)   = 0.01879d0 !Bratmixave

        ABtransmission(1,1) = 0.4841d0 !AdMCave
        ABtransmission(2,1) = 0.3632d0 !AdMCdev
        ABtransmission(1,2) = 0.4834d0 !AdLPave
        ABtransmission(1,3) = 0.48537d0 !BrMCave
        ABtransmission(1,4) = 0.48446d0 !BrLPave
        ABtransmission(1,5) = 0.23065d0 !Bratmixave
      elseif( Adamscase<2.3d0+eps .and. Adamscase>2.3d0-eps ) then
        ABreflection(1,1)   = 0.0856d0 !AdMCave
        ABreflection(2,1)   = 0.0697d0 !AdMCdev
        ABreflection(1,2)   = 0.0591d0 !AdLPave
        ABreflection(1,3)   = 0.08549d0 !BrMCave
        ABreflection(1,4)   = 0.05864d0 !BrLPave
        ABreflection(1,5)   = 0.01963d0 !Bratmixave

        ABtransmission(1,1) = 0.0016d0 !AdMCave
        ABtransmission(2,1) = 0.0121d0 !AdMCdev
        ABtransmission(1,2) = 0.0015d0 !AdLPave
        ABtransmission(1,3) = 0.00166d0 !BrMCave
        ABtransmission(1,4) = 0.00154d0 !BrLPave
        ABtransmission(1,5) = 0.00001d0 !Bratmixave
      elseif( Adamscase<3.1d0+eps .and. Adamscase>3.1d0-eps ) then
        ABreflection(1,1)   = 0.0480d0 !AdMCave
        ABreflection(2,1)   = 0.0935d0 !AdMCdev
        ABreflection(1,2)   = 0.0473d0 !AdLPave
        ABreflection(1,3)   = 0.04812d0 !BrMCave
        ABreflection(1,4)   = 0.04702d0 !BrLPave
        ABreflection(1,5)   = 0.07455d0 !Bratmixave

        ABtransmission(1,1) = 0.9341d0 !AdMCave
        ABtransmission(2,1) = 0.1339d0 !AdMCdev
        ABtransmission(1,2) = 0.9344d0 !AdLPave
        ABtransmission(1,3) = 0.93394d0 !BrMCave
        ABtransmission(1,4) = 0.93468d0 !BrLPave
        ABtransmission(1,5) = 0.90602d0 !Bratmixave
      elseif( Adamscase<3.2d0+eps .and. Adamscase>3.2d0-eps ) then
        ABreflection(1,1)   = 0.2563d0 !AdMCave
        ABreflection(2,1)   = 0.1546d0 !AdMCdev
        ABreflection(1,2)   = 0.2178d0 !AdLPave
        ABreflection(1,3)   = 0.25419d0 !BrMCave
        ABreflection(1,4)   = 0.21677d0 !BrLPave
        ABreflection(1,5)   = 0.35288d0 !Bratmixave

        ABtransmission(1,1) = 0.5985d0 !AdMCave
        ABtransmission(2,1) = 0.2832d0 !AdMCdev
        ABtransmission(1,2) = 0.6267d0 !AdLPave
        ABtransmission(1,3) = 0.60114d0 !BrMCave
        ABtransmission(1,4) = 0.62755d0 !BrLPave
        ABtransmission(1,5) = 0.47477d0 !Bratmixave
      elseif( Adamscase<3.3d0+eps .and. Adamscase>3.3d0-eps ) then
        ABreflection(1,1)   = 0.4785d0 !AdMCave
        ABreflection(2,1)   = 0.0038d0 !AdMCdev
        ABreflection(1,2)   = 0.3707d0 !AdLPave
        ABreflection(1,3)   = 0.47746d0 !BrMCave
        ABreflection(1,4)   = 0.36952d0 !BrLPave
        ABreflection(1,5)   = 0.47819d0 !Bratmixave

        ABtransmission(1,1) = 0.0159d0 !AdMCave
        ABtransmission(2,1) = 0.0314d0 !AdMCdev
        ABtransmission(1,2) = 0.0237d0 !AdLPave
        ABtransmission(1,3) = 0.01609d0 !BrMCave
        ABtransmission(1,4) = 0.02377d0 !BrLPave
        ABtransmission(1,5) = 0.00386d0 !Bratmixave
      elseif( Adamscase<4.1d0+eps .and. Adamscase>4.1d0-eps ) then
        ABreflection(1,1)   = 0.0434d0 !AdMCave
        ABreflection(2,1)   = 0.1267d0 !AdMCdev
        ABreflection(1,2)   = 0.0432d0 !AdLPave
        ABreflection(1,3)   = 0.04304d0 !BrMCave
        ABreflection(1,4)   = 0.04302d0 !BrLPave
        ABreflection(1,5)   = 0.07544d0 !Bratmixave

        ABtransmission(1,1) = 0.9388d0 !AdMCave
        ABtransmission(2,1) = 0.1209d0 !AdMCdev
        ABtransmission(1,2) = 0.9390d0 !AdLPave
        ABtransmission(1,3) = 0.93940d0 !BrMCave
        ABtransmission(1,4) = 0.93931d0 !BrLPave
        ABtransmission(1,5) = 0.90690d0 !Bratmixave
      elseif( Adamscase<4.2d0+eps .and. Adamscase>4.2d0-eps ) then
        ABreflection(1,1)   = 0.1224d0 !AdMCave
        ABreflection(2,1)   = 0.2726d0 !AdMCdev
        ABreflection(1,2)   = 0.1068d0 !AdLPave
        ABreflection(1,3)   = 0.12120d0 !BrMCave
        ABreflection(1,4)   = 0.10681d0 !BrLPave
        ABreflection(1,5)   = 0.36040d0 !Bratmixave

        ABtransmission(1,1) = 0.7233d0 !AdMCave
        ABtransmission(2,1) = 0.2306d0 !AdMCdev
        ABtransmission(1,2) = 0.7385d0 !AdLPave
        ABtransmission(1,3) = 0.72656d0 !BrMCave
        ABtransmission(1,4) = 0.74079d0 !BrLPave
        ABtransmission(1,5) = 0.48096d0 !Bratmixave
      elseif( Adamscase<4.3d0+eps .and. Adamscase>4.3d0-eps ) then
        ABreflection(1,1)   = 0.2369d0 !AdMCave
        ABreflection(2,1)   = 0.2860d0 !AdMCdev
        ABreflection(1,2)   = 0.1799d0 !AdLPave
        ABreflection(1,3)   = 0.23723d0 !BrMCave
        ABreflection(1,4)   = 0.18049d0 !BrLPave
        ABreflection(1,5)   = 0.49564d0 !Bratmixave

        ABtransmission(1,1) = 0.0981d0 !AdMCave
        ABtransmission(2,1) = 0.0887d0 !AdMCdev
        ABtransmission(1,2) = 0.1278d0 !AdLPave
        ABtransmission(1,3) = 0.09843d0 !BrMCave
        ABtransmission(1,4) = 0.12840d0 !BrLPave
        ABtransmission(1,5) = 0.00474d0 !Bratmixave
      elseif( Adamscase<5.1d0+eps .and. Adamscase>5.1d0-eps ) then
        ABreflection(1,1)   = 0.0089d0 !AdMCave
        ABreflection(2,1)   = 0.0030d0 !AdMCdev
        ABreflection(1,2)   = 0.0089d0 !AdLPave
        ABreflection(1,3)   = 0.00885d0 !BrMCave
        ABreflection(1,4)   = 0.00880d0 !BrLPave
        ABreflection(1,5)   = 0.00654d0 !Bratmixave

        ABtransmission(1,1) = 0.9140d0 !AdMCave
        ABtransmission(2,1) = 0.2217d0 !AdMCdev
        ABtransmission(1,2) = 0.9140d0 !AdLPave
        ABtransmission(1,3) = 0.91465d0 !BrMCave
        ABtransmission(1,4) = 0.91419d0 !BrLPave
        ABtransmission(1,5) = 0.83900d0 !Bratmixave
      elseif( Adamscase<5.2d0+eps .and. Adamscase>5.2d0-eps ) then
        ABreflection(1,1)   = 0.0744d0 !AdMCave
        ABreflection(2,1)   = 0.0278d0 !AdMCdev
        ABreflection(1,2)   = 0.0717d0 !AdLPave
        ABreflection(1,3)   = 0.07337d0 !BrMCave
        ABreflection(1,4)   = 0.07063d0 !BrLPave
        ABreflection(1,5)   = 0.01879d0 !Bratmixave

        ABtransmission(1,1) = 0.7588d0 !AdMCave
        ABtransmission(2,1) = 0.3325d0 !AdMCdev
        ABtransmission(1,2) = 0.7581d0 !AdLPave
        ABtransmission(1,3) = 0.75974d0 !BrMCave
        ABtransmission(1,4) = 0.75920d0 !BrLPave
        ABtransmission(1,5) = 0.23065d0 !Bratmixave
      elseif( Adamscase<5.3d0+eps .and. Adamscase>5.3d0-eps ) then
        ABreflection(1,1)   = 0.2897d0 !AdMCave
        ABreflection(2,1)   = 0.1631d0 !AdMCdev
        ABreflection(1,2)   = 0.2193d0 !AdLPave
        ABreflection(1,3)   = 0.28763d0 !BrMCave
        ABreflection(1,4)   = 0.21829d0 !BrLPave
        ABreflection(1,5)   = 0.01963d0 !Bratmixave

        ABtransmission(1,1) = 0.1960d0 !AdMCave
        ABtransmission(2,1) = 0.2551d0 !AdMCdev
        ABtransmission(1,2) = 0.1787d0 !AdLPave
        ABtransmission(1,3) = 0.19553d0 !BrMCave
        ABtransmission(1,4) = 0.17938d0 !BrLPave
        ABtransmission(1,5) = 0.00001d0 !Bratmixave
      elseif( Adamscase<6.1d0+eps .and. Adamscase>6.1d0-eps ) then
        ABreflection(1,1)   = 0.0426d0 !AdMCave
        ABreflection(2,1)   = 0.0985d0 !AdMCdev
        ABreflection(1,2)   = 0.0426d0 !AdLPave
        ABreflection(1,3)   = 0.04246d0 !BrMCave
        ABreflection(1,4)   = 0.04236d0 !BrLPave
        ABreflection(1,5)   = 0.07455d0 !Bratmixave

        ABtransmission(1,1) = 0.9398d0 !AdMCave
        ABtransmission(2,1) = 0.1446d0 !AdMCdev
        ABtransmission(1,2) = 0.9397d0 !AdLPave
        ABtransmission(1,3) = 0.94005d0 !BrMCave
        ABtransmission(1,4) = 0.93985d0 !BrLPave
        ABtransmission(1,5) = 0.90602d0 !Bratmixave
      elseif( Adamscase<6.2d0+eps .and. Adamscase>6.2d0-eps ) then
        ABreflection(1,1)   = 0.1440d0 !AdMCave
        ABreflection(2,1)   = 0.1452d0 !AdMCdev
        ABreflection(1,2)   = 0.1255d0 !AdLPave
        ABreflection(1,3)   = 0.14248d0 !BrMCave
        ABreflection(1,4)   = 0.12426d0 !BrLPave
        ABreflection(1,5)   = 0.35288d0 !Bratmixave

        ABtransmission(1,1) = 0.7666d0 !AdMCave
        ABtransmission(2,1) = 0.3011d0 !AdMCdev
        ABtransmission(1,2) = 0.7733d0 !AdLPave
        ABtransmission(1,3) = 0.76855d0 !BrMCave
        ABtransmission(1,4) = 0.77434d0 !BrLPave
        ABtransmission(1,5) = 0.47477d0 !Bratmixave
      elseif( Adamscase<6.3d0+eps .and. Adamscase>6.3d0-eps ) then
        ABreflection(1,1)   = 0.4344d0 !AdMCave
        ABreflection(2,1)   = 0.0572d0 !AdMCdev
        ABreflection(1,2)   = 0.2910d0 !AdLPave
        ABreflection(1,3)   = 0.43319d0 !BrMCave
        ABreflection(1,4)   = 0.28962d0 !BrLPave
        ABreflection(1,5)   = 0.47819d0 !Bratmixave

        ABtransmission(1,1) = 0.1861d0 !AdMCave
        ABtransmission(2,1) = 0.2134d0 !AdMCdev
        ABtransmission(1,2) = 0.1945d0 !AdLPave
        ABtransmission(1,3) = 0.18690d0 !BrMCave
        ABtransmission(1,4) = 0.19498d0 !BrLPave
        ABtransmission(1,5) = 0.00386d0 !Bratmixave
      elseif( Adamscase<7.1d0+eps .and. Adamscase>7.1d0-eps ) then
        ABreflection(1,1)   = 0.0763d0 !AdMCave
        ABreflection(2,1)   = 0.0751d0 !AdMCdev
        ABreflection(1,2)   = 0.0758d0 !AdLPave
        ABreflection(1,3)   = 0.07494d0 !BrMCave
        ABreflection(1,4)   = 0.07491d0 !BrLPave
        ABreflection(1,5)   = 0.08343d0 !Bratmixave

        ABtransmission(1,1) = 0.9218d0 !AdMCave
        ABtransmission(2,1) = 0.0732d0 !AdMCdev
        ABtransmission(1,2) = 0.9223d0 !AdLPave
        ABtransmission(1,3) = 0.92329d0 !BrMCave
        ABtransmission(1,4) = 0.92312d0 !BrLPave
        ABtransmission(1,5) = 0.91479d0 !Bratmixave
      elseif( Adamscase<7.2d0+eps .and. Adamscase>7.2d0-eps ) then
        ABreflection(1,1)   = 0.3210d0 !AdMCave
        ABreflection(2,1)   = 0.2865d0 !AdMCdev
        ABreflection(1,2)   = 0.3157d0 !AdLPave
        ABreflection(1,3)   = 0.32096d0 !BrMCave
        ABreflection(1,4)   = 0.31515d0 !BrLPave
        ABreflection(1,5)   = 0.43624d0 !Bratmixave

        ABtransmission(1,1) = 0.6599d0 !AdMCave
        ABtransmission(2,1) = 0.2688d0 !AdMCdev
        ABtransmission(1,2) = 0.6652d0 !AdLPave
        ABtransmission(1,3) = 0.66029d0 !BrMCave
        ABtransmission(1,4) = 0.66595d0 !BrLPave
        ABtransmission(1,5) = 0.54446d0 !Bratmixave
      elseif( Adamscase<7.3d0+eps .and. Adamscase>7.3d0-eps ) then
        ABreflection(1,1)   = 0.6916d0 !AdMCave
        ABreflection(2,1)   = 0.2615d0 !AdMCdev
        ABreflection(1,2)   = 0.6070d0 !AdLPave
        ABreflection(1,3)   = 0.69109d0 !BrMCave
        ABreflection(1,4)   = 0.60759d0 !BrLPave
        ABreflection(1,5)   = 0.78629d0 !Bratmixave

        ABtransmission(1,1) = 0.1615d0 !AdMCave
        ABtransmission(2,1) = 0.1740d0 !AdMCdev
        ABtransmission(1,2) = 0.2391d0 !AdLPave
        ABtransmission(1,3) = 0.16350d0 !BrMCave
        ABtransmission(1,4) = 0.24038d0 !BrLPave
        ABtransmission(1,5) = 0.06677d0 !Bratmixave
      elseif( Adamscase<8.1d0+eps .and. Adamscase>8.1d0-eps ) then
        ABreflection(1,1)   = 0.0010d0 !AdMCave
        ABreflection(2,1)   = 0.0010d0 !AdMCdev
        ABreflection(1,2)   = 0.0010d0 !AdLPave
        ABreflection(1,3)   = 0.00098d0 !BrMCave
        ABreflection(1,4)   = 0.00098d0 !BrLPave
        ABreflection(1,5)   = 0.00070d0 !Bratmixave

        ABtransmission(1,1) = 0.8509d0 !AdMCave
        ABtransmission(2,1) = 0.1464d0 !AdMCdev
        ABtransmission(1,2) = 0.8503d0 !AdLPave
        ABtransmission(1,3) = 0.85192d0 !BrMCave
        ABtransmission(1,4) = 0.85188d0 !BrLPave
        ABtransmission(1,5) = 0.83326d0 !Bratmixave
      elseif( Adamscase<8.2d0+eps .and. Adamscase>8.2d0-eps ) then
        ABreflection(1,1)   = 0.0088d0 !AdMCave
        ABreflection(2,1)   = 0.0091d0 !AdMCdev
        ABreflection(1,2)   = 0.0085d0 !AdLPave
        ABreflection(1,3)   = 0.00877d0 !BrMCave
        ABreflection(1,4)   = 0.00840d0 !BrLPave
        ABreflection(1,5)   = 0.00196d0 !Bratmixave

        ABtransmission(1,1) = 0.4818d0 !AdMCave
        ABtransmission(2,1) = 0.4361d0 !AdMCdev
        ABtransmission(1,2) = 0.4826d0 !AdLPave
        ABtransmission(1,3) = 0.48297d0 !BrMCave
        ABtransmission(1,4) = 0.48293d0 !BrLPave
        ABtransmission(1,5) = 0.22054d0 !Bratmixave
      elseif( Adamscase<8.3d0+eps .and. Adamscase>8.3d0-eps ) then
        ABreflection(1,1)   = 0.0369d0 !AdMCave
        ABreflection(2,1)   = 0.0500d0 !AdMCdev
        ABreflection(1,2)   = 0.0243d0 !AdLPave
        ABreflection(1,3)   = 0.03651d0 !BrMCave
        ABreflection(1,4)   = 0.02402d0 !BrLPave
        ABreflection(1,5)   = 0.00204d0 !Bratmixave

        ABtransmission(1,1) = 0.0766d0 !AdMCave
        ABtransmission(2,1) = 0.2252d0 !AdMCdev
        ABtransmission(1,2) = 0.0755d0 !AdLPave
        ABtransmission(1,3) = 0.07678d0 !BrMCave
        ABtransmission(1,4) = 0.07567d0 !BrLPave
        ABtransmission(1,5) = 0.00001d0 !Bratmixave
      elseif( Adamscase<9.1d0+eps .and. Adamscase>9.1d0-eps ) then
        ABreflection(1,1)   = 0.0670d0 !AdMCave
        ABreflection(2,1)   = 0.0646d0 !AdMCdev
        ABreflection(1,2)   = 0.0669d0 !AdLPave
        ABreflection(1,3)   = 0.06618d0 !BrMCave
        ABreflection(1,4)   = 0.06608d0 !BrLPave
        ABreflection(1,5)   = 0.07455d0 !Bratmixave

        ABtransmission(1,1) = 0.9136d0 !AdMCave
        ABtransmission(2,1) = 0.0834d0 !AdMCdev
        ABtransmission(1,2) = 0.9137d0 !AdLPave
        ABtransmission(1,3) = 0.91468d0 !BrMCave
        ABtransmission(1,4) = 0.91458d0 !BrLPave
        ABtransmission(1,5) = 0.90602d0 !Bratmixave
      elseif( Adamscase<9.2d0+eps .and. Adamscase>9.2d0-eps ) then
        ABreflection(1,1)   = 0.2435d0 !AdMCave
        ABreflection(2,1)   = 0.1991d0 !AdMCdev
        ABreflection(1,2)   = 0.2381d0 !AdLPave
        ABreflection(1,3)   = 0.24286d0 !BrMCave
        ABreflection(1,4)   = 0.23727d0 !BrLPave
        ABreflection(1,5)   = 0.35288d0 !Bratmixave

        ABtransmission(1,1) = 0.6045d0 !AdMCave
        ABtransmission(2,1) = 0.3346d0 !AdMCdev
        ABtransmission(1,2) = 0.6086d0 !AdLPave
        ABtransmission(1,3) = 0.60498d0 !BrMCave
        ABtransmission(1,4) = 0.60911d0 !BrLPave
        ABtransmission(1,5) = 0.47477d0 !Bratmixave
      elseif( Adamscase<9.3d0+eps .and. Adamscase>9.3d0-eps ) then
        ABreflection(1,1)   = 0.4466d0 !AdMCave
        ABreflection(2,1)   = 0.0923d0 !AdMCdev
        ABreflection(1,2)   = 0.3272d0 !AdLPave
        ABreflection(1,3)   = 0.44516d0 !BrMCave
        ABreflection(1,4)   = 0.32612d0 !BrLPave
        ABreflection(1,5)   = 0.47819d0 !Bratmixave

        ABtransmission(1,1) = 0.1037d0 !AdMCave
        ABtransmission(2,1) = 0.2290d0 !AdMCdev
        ABtransmission(1,2) = 0.1195d0 !AdLPave
        ABtransmission(1,3) = 0.10457d0 !BrMCave
        ABtransmission(1,4) = 0.11967d0 !BrLPave
        ABtransmission(1,5) = 0.00386d0 !Bratmixave
      endif
    endif

  end subroutine Acase_load

end module Loadcase

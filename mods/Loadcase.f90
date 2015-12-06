module Loadcase
  implicit none

CONTAINS
  ! print statements in this module use # 100-199

  subroutine read_test_inputstoc
  use rngvars, only: rngseed
  use genRealzvars,         only: Adamscase, sig, scatrat, lam, s, numRealz, pltgenrealznumof, &
                                  pltgenrealz, pltgenrealzwhich, GBsigave, GBsigvar, GBscatrat, &
                                  GBlamc, GBs, flCorrMarkov, flCorrRealz, chgeomtype
  use KLvars,               only: KLvarcalc, KLvarkept_tol, pltEigfwhich, pltxiBinswhich, &
                                  pltCowhich, pltxiBinsnumof, pltEigfnumof, pltConumof, binNumof,&
                                  numEigs, numSlice, levsrefEig, Corrnumpoints, binSmallBound, &
                                  binLargeBound, pltxiBins, pltxiBinsgauss, pltEigf, pltCo, &
                                  Corropts, KLrnumpoints, pltKLrealz, pltKLrealznumof, pltKLrealzwhich, &
                                  flmeanadjust, meanadjust_tol, &
                                  Gaussrandtype, flCorrKL, numrefinesameiter, flGaussdiffrand, &
                                  chGausstype, chLNmode, flLNxscheck, numLNxspts, numLNxsbins, &
                                  chLNxschecktype, chLNxsplottype
  use MCvars,               only: trprofile_binnum, binplot, numParts, trannprt, rodOrplanar, sourceType, &
                                  pltflux, flnegxs, LPamnumParts, fluxnumcells, pltmatflux, &
                                  pltfluxtype, flCorrMC, chTrantype
  character(7) :: pltallopt                         !Plot all same opt
  character(3) :: allL1s,allcenters !all basic L1/center based functionals?
  character(6) :: chosennorm        !norm to convert from text to number

  real(8)       :: dumreal,s2
  character(4)  :: setflags(3)
  character(20) :: dumchar !use this to "skip" a line
  integer         :: i, ones, zeros
  real(8)         :: eps = 0.000001d0
  logical         :: flstopstatus = .false., flsleep = .false.

  open(unit=2,file="inputstoc.txt")

  read(2,*) rngseed
  read(2,*) chgeomtype

  !--- Geometry - num of realz ---!
  read(2,*) dumchar
  read(2,*) numRealz,trannprt

  !--- Geometry - Gauss or Gauss-based type problem ---!
  read(2,*) dumchar
  read(2,*) chGausstype
  read(2,*) chLNmode,setflags(1)
  if(setflags(1)=='same') flGaussdiffrand = .false.
  read(2,*) GBsigave,GBsigvar
  read(2,*) GBscatrat
  read(2,*) GBlamc
  read(2,*) GBs

  !--- Geometry - 'Markov' type problem ---!
  read(2,*) dumchar
  read(2,*) Adamscase
  read(2,*) sig(1),sig(2)
  read(2,*) scatrat(1),scatrat(2)
  read(2,*) lam(1),lam(2)
  read(2,*) s
  read(2,*) setflags(1),setflags(2),setflags(3)
  if(setflags(1)=='yes') flCorrMarkov=.true.
  if(setflags(2)=='yes') flCorrKL    =.true.
  if(setflags(3)=='yes') flCorrRealz =.true.

  !--- Large KL Options ---!
  read(2,*) dumchar
  read(2,*) numEigs
  read(2,*) binNumof

  !--- Large MCtrans Options ---!
  read(2,*) dumchar
  read(2,*) chTrantype
  read(2,*) numParts
  read(2,*) LPamnumParts
  read(2,*) setflags(1)
  if(setflags(1)=='yes') flCorrMC    =.true.

  !--- Lesser KL Options ---!
  read(2,*) dumchar
  read(2,*) KLrnumpoints
  read(2,*) levsrefEig
  read(2,*) numrefinesameiter
  read(2,*) binSmallBound,binLargeBound
  read(2,*) KLvarcalc,KLvarkept_tol
  read(2,*) numSlice
  read(2,*) Gaussrandtype

  !--- Lesser MCtrans Options ---!
  read(2,*) dumchar 
  read(2,*) rodOrplanar
  read(2,*) sourceType
  read(2,*) setflags(1)
  if(setflags(1)=='yes') flnegxs  =.true.
  read(2,*) setflags(1),meanadjust_tol
  if(setflags(1)=='yes') flmeanadjust=.true.
  read(2,*) setflags(1),chLNxschecktype,numLNxspts,numLNxsbins,chLNxsplottype
  if(setflags(1)=='yes') flLNxscheck=.true.


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
  allocate(pltKLrealzwhich(2,pltKLrealznumof))
  do i=1,pltKLrealznumof
    read(2,*) pltKLrealzwhich(1,i),pltKLrealzwhich(2,i)
  enddo

  read(2,*) dumchar    !Plotting Variace (Co)
  read(2,*) pltCo(1),pltCo(2),pltCo(3),pltCo(4)
  read(2,*) pltConumof
  allocate(pltCowhich(2,pltConumof))
  do i=1,pltConumof
    read(2,*) pltCowhich(1,i),pltCowhich(2,i)
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




  !begin tests of valid input
  print *,"  "

                      !Tests for geometry vs transport type
  if(chgeomtype=='contin' .and. .not.(chTrantype=='GaussKL' .or. chTrantype=='None')) then
    print *,"--User attempting to use invalid transport type for continuous geometry"
    flstopstatus = .true.
  endif
  if(chgeomtype=='binary' .and. .not.(chTrantype=='radMC' .or. chTrantype=='radWood' .or.&
     chTrantype=='KLWood'.or. chTrantype=='LPMC'  .or. chTrantype=='atmixMC' .or. chTrantype=='None')) then
    print *,"--User attempting to use invalid transport type for binary material geometry"
    flstopstatus = .true.
  endif

  do i=1,pltEigfnumof    !Test Eigenfunction plotting order of Eigs
    if( pltEigfwhich(i)>numEigs .AND. pltEigf(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot Eigenfunction of higher order than Eigenvalues calculated"
      flstopstatus = .true.
    endif
  enddo

  do i=1,pltxiBinsnumof  !Test xiBins plotting order of Eigs and num of bins
    if( pltxiBinswhich(1,i)>numEigs .AND. pltxiBins(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot xiBins for Eigenvalue of higher order than calculated"
      flstopstatus = .true.
    endif
    if( pltxiBinswhich(2,i)>numRealz .AND. pltxiBins(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot xibins for more realizations than generated"
      flstopstatus = .true.
    endif
  enddo

  do i=1,pltKLrealznumof
    if( pltKLrealzwhich(1,i)>numRealz .AND. pltKLrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot more reconstructed realz than reconstructed"
      flstopstatus = .true.
    endif
    if( pltKLrealzwhich(2,i)>numEigs .AND. pltKLrealz(1) .NE. 'noplot' ) then
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

                              !Test radtransMC print frequency
  if( (trannprt>numRealz .or. mod(numRealz,trannprt)/=0) .and. chTrantype=='radMC' ) then 
    print *,"--Print to screen frequency for MCtran must be factor of number of realizations"
    flstopstatus = .true.
  endif

                              !Test plotting flux
  if( pltflux(1)/='noplot' .and. .not.chTrantype=='radMC' .and. .not.chTrantype=='radWood' &
      .and. .not.chTrantype=='KLWood' .and. .not.chTrantype=='LPMC' .and. &
      .not.chTrantype=='atmixMC' .and. .not.chTrantype=='GaussKL') then
    print *,"--User attempting to plot flux when no transport calculations are made"
    flstopstatus = .true.
  endif
  if( sourceType/='left' .AND. sourceType/='intern' ) then
    print *,"--User attempting to run invalid source type.  Please put either 'left' or 'intern'"
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
    if( pltCowhich(1,i)>numEigs .AND. pltCo(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot CoEff values for eigenvalues not calculated"
      flstopstatus = .true.
    endif
    if( pltCowhich(2,i)/=1 .AND. pltCowhich(2,i)/=2 .AND. pltCo(1) .NE. 'noplot' ) then
      print *,"--User input for 'CoEffExp vs CoEffAct' not valid.  Enter a '1' or a '2'"
      flstopstatus = .true.
    endif
  enddo

  if( flCorrRealz .and. (.not.flCorrMarkov .or. .not.flCorrKL)) then !Test for correlation deficiency
    print *,"--User trying to correlate KLres w/ others w/o corr Markov or KL, both set to true"
    flCorrMarkov=.true.
    flCorrKL    =.true.
    flsleep = .true.
  endif

  if( Gaussrandtype=='BM' .and. flCorrKL ) then  !Test for correlated GaussKL sampling method
    print *,"--User trying to correlate KL realz w/ Box-Muller sampling, switched to inverse sampling"
    Gaussrandtype='inv'
    flsleep = .true.
  endif

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

  end subroutine read_test_inputstoc




  subroutine global_allocate
  !This subroutine allocates and initializes all global variables
  use rngvars, only: rngappnum, rngseed
  use genRealzvars, only: lam, P, s, numRealz, numPath, sumPath, sqrPath, largesti, &
                          totLength, lamc, sig, sigave, sigscatave, sigabsave, scatrat, &
                          numPosRealz, numNegRealz, numRealz, sigvar, atmixscatrat, &
                          scatvar, absvar, atmixsig, chgeomtype, GBs, GBlamc, GBscatrat, &
                          GBsigvar, GBsigave
  use KLvars, only: KLrnumpoints, numEigs, pltKLrealznumof, chGausstype, Corropts, &
                    KLrxisig, numSlice, gam, alpha, Ak, Eig, chLNmode, pltCo, &
                    xi, KLrxivals, KLrxivalss, pltKLrealzarray, flGaussdiffrand, pltKLrealz
  use MCvars, only: fluxfaces, numParts, stocMC_reflection, stocMC_transmission, &
                    stocMC_absorption, LPamnumParts, stocMC_fluxall, chTrantype, &
                    stocMC_fluxmat1, stocMC_fluxmat2, pltflux, pltmatflux, areapnsamp, &
                    fluxnumcells, flfluxplot, LPamMCsums, transmit, reflect, absorb, &
                    numpnSamp, radtrans_int, Wood_rej
  use mcnp_random, only: RN_init_problem
  use utilities, only: exponentialfit
  integer :: i

  !initialize rngvars
  rngappnum  = 0
  call RN_init_problem( 1, rngseed, int(0,8), int(0,8), 0)


  !allocate and initialize genRealzvars
  numPosRealz= 0
  numNegRealz= 0
  if(chgeomtype=='contin') then  !Gauss-based input
    sigave       = GBsigave
    sigvar       = GBsigvar
    lamc         = GBlamc
    s            = GBs
    if(chGausstype=='LogN') then
      sigave  = log(GBsigave**2/sqrt(sigvar+GBsigave**2))
      sigvar  = log(sigvar/GBsigave**2+1.0d0)
      sig(1)  = log(sig(1)**2/sqrt(sigvar+sig(1) **2))
      sig(2)  = log(sig(2)**2/sqrt(sigvar+sig(2) **2))
      if(chLNmode=='fitlamc') lamc = exponentialfit(s,1d0+sigvar/sigave,lamc)
    endif
    sigscatave = sigave *      GBscatrat
    sigabsave  = sigave * (1d0-GBscatrat)
    scatvar    = sigvar *      GBscatrat
    absvar     = sigvar * (1d0-GBscatrat)
  elseif(chgeomtype=='binary') then
    numPath    = 0  !setup Markov material tallies
    sumPath    = 0d0
    sqrPath    = 0d0
    largesti   = 0
    totLength  = 0d0
    P(1)       = lam(1)/(lam(1)+lam(2)) !calc probabilities
    P(2)       = lam(2)/(lam(1)+lam(2))
    lamc       = (lam(1)*lam(2))/(lam(1)+lam(2))
    sigave     = P(1)*                 sig(1) + P(2)*                 sig(2)
    sigscatave = P(1)*     scatrat(1) *sig(1) + P(2)*     scatrat(2) *sig(2)
    sigabsave  = P(1)*(1d0-scatrat(1))*sig(1) + P(2)*(1d0-scatrat(2))*sig(2)
    sigvar     = P(1)*P(2) * (sig(1)                  - sig(2)                ) **2
    scatvar    = P(1)*P(2) * (sig(1)*     scatrat(1)  - sig(2)*     scatrat(2)) **2
    absvar     = P(1)*P(2) * (sig(1)*(1d0-scatrat(1)) - sig(2)*(1d0-scatrat(2)))**2
    if(chTrantype=='atmixMC') then
      atmixsig     =   P(1)*sig(1)            + P(2)*sig(2)
      atmixscatrat = ( P(1)*sig(1)*scatrat(1) + P(2)*sig(2)*scatrat(2) ) / atmixsig
    endif
  endif


  !allocate  KLresearch variables
  if(chTrantype=='KLWood' .or. chTrantype=='GaussKL' .or. &
     Corropts(1).ne.'noplot' .or. pltCo(1).ne.'noplot' .or. pltKLrealz(1).ne.'noplot') then
    allocate(gam(numEigs))
    allocate(alpha(numEigs))
    allocate(Ak(numEigs))
    allocate(Eig(numEigs))
    allocate(xi(numRealz,numEigs))
  endif


  !allocate and initialize KLconstruction variables
  if(chTrantype=='KLWood' .or. chTrantype=='GaussKL' .or. pltKLrealz(1).ne.'noplot') then
    allocate(KLrxivals(numRealz,numEigs))
    if(chTrantype=='GaussKL' .and. flGaussdiffrand) allocate(KLrxivalss(numRealz,numEigs))
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
    allocate(stocMC_reflection(2))   !global MC variables for each method
    allocate(stocMC_transmission(2)) !rank 2 holds 1=average, 2=deviation
    allocate(stocMC_absorption(2))
    stocMC_reflection   = 0.0d0
    stocMC_transmission = 0.0d0
    stocMC_absorption   = 0.0d0
    if(.not.allocated(transmit)) allocate(transmit(numRealz)) !leakage tallies
    if(.not.allocated(reflect))  allocate(reflect(numRealz))
    if(.not.allocated(absorb))   allocate(absorb(numRealz))
    transmit     = 0.0d0
    reflect      = 0.0d0
    absorb       = 0.0d0

    flfluxplot = .false.  !flux variable allocations
    if(.not.(pltflux(1)=='noplot' .and. pltmatflux=='noplot')) flfluxplot = .true.
    if(flfluxplot) then   !alloc and init flux cells
      allocate(fluxfaces(fluxnumcells+1))
      fluxfaces = 0.0d0
      do i=1,fluxnumcells+1
        fluxfaces(i) = (s/fluxnumcells) * (i-1)
      enddo
    endif
    if(.not.pltflux(1)=='noplot' .or. &!mat irrespective flux allocations
      (flfluxplot .and. (chTrantype=='LPMC' .or. chTrantype=='atmixMC')) ) then
      allocate(stocMC_fluxall(fluxnumcells,2))
      stocMC_fluxall = 0.0d0
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
  endif

  end subroutine global_allocate




  subroutine clearreports
  call system("rm texts/*.out")
  end subroutine clearreports



  subroutine finalreport
  call system("cat texts/Acase.out texts/Woodnegstats.out texts/MCleakage.out > texts/finalreport.out")
  call system("cat texts/finalreport.out")
  end subroutine finalreport



  subroutine Acase_load
  use genRealzvars, only:   Adamscase, sig, lam, scatrat, s  
  use MCvars, only: rodOrplanar, ABreflection, ABtransmission
  real(8) :: eps = 0.0001d0 !local var
    ! Loading specials cases from Adams and Pomraning Paper, 27 total cases
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
        s=0.1d0
      elseif( int((Adamscase*10+0.005d0-(int(Adamscase))*10))==2 ) then
        s=1.0d0
      elseif( int((Adamscase*10+0.005d0-(int(Adamscase))*10))==3 ) then
        s=10.0d0
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



  subroutine Acase_print
  use genRealzvars, only: Adamscase, sig, lam, scatrat, s, lamc

  100 format("  Values loaded for special Adamscase :",f5.1)
  101 format("    sig(1)    :",f10.6,"   sig(2)    :",f10.6)
  102 format("    scatrat(1):",f6.2,"       scatrat(2):",f6.2)
  103 format("    lam(1)    :",f6.2,"       lam(2)    :",f6.2)
  104 format("    s         :",f6.2,"       lamc      :",f6.2)

  open(unit=100,file="Acase.out")
  write(100,*)
  write(100,100) Adamscase
  write(100,101) sig(1),sig(2)
  write(100,102) scatrat(1),scatrat(2)
  write(100,103) lam(1),lam(2)
  write(100,104) s,lamc
  write(100,*)

  close(unit=100)
  call system("mv Acase.out texts")

  end subroutine Acase_print



end module Loadcase

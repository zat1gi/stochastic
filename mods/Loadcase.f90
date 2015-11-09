module Loadcase
  implicit none

CONTAINS
  ! print statements in this module use # 100-199



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





  subroutine readinputstoc
  use rngvars, only: rngseed
  use genRealzvars,         only: Adamscase, sig, scatrat, lam, s, numRealz, pltgenrealznumof, &
                                  pltgenrealz, pltgenrealzwhich, GBsigave, GBsigvar, GBscatrat, &
                                  GBlamc, GBs, flCorrMarkov, flCorrRealz, flGBgeom
  use genSampvars, only: specialprob, nummat, param1, param2, param1_mean, param1_uncert, &
                         param2_mean, param2_uncert
  use KLvars,               only: KLvarcalc, KLvarkept_tol, pltEigfwhich, pltxiBinswhich, &
                                  pltCowhich, pltxiBinsnumof, pltEigfnumof, pltConumof, binNumof,&
                                  numEigs, numSlice, levsrefEig, Corrnumpoints, binSmallBound, &
                                  binLargeBound, pltxiBins, pltxiBinsgauss, pltEigf, pltCo, &
                                  Corropts, KLrnumpoints, KLrnumRealz, KLrprintat, pltKLrrealz, &
                                  pltKLrrealznumof, pltKLrrealzwhich, pltKLrrealzPointorXi, &
                                  KLres, KLrec, KLnoise, flmatbasedxs, flmeanadjust, meanadjust_tol, &
                                  Gaussrandtype, flCorrKL, numrefinesameiter, flglGaussdiffrand, &
                                  flglLN, chLNmode, flLNxscheck, numLNxspts, numLNxsbins
  use MCvars,               only: trprofile_binnum, radMCbinplot, radWoodbinplot, KLWoodbinplot, &
                                  GaussKLbinplot, numParts, trannprt, rodOrplanar, sourceType, &
                                  pltflux, flnegxs, fldistneg, radMC, radWood, WAMC, GaussKL, &
                                  KLWood, LPMC, atmixMC, LPamnumParts, fluxnumcells, pltmatflux, &
                                  pltfluxtype, refsigMode, userrefsig, wgtmax, wgtmin, wgtmaxmin, &
                                  negwgtbinnum, nwvalsperbin, probtype, flCorrMC
  use MLMCvars,             only: detMLMC, MLMC_TOL, numcellsLevel0, nextLevelFactor, MLMC_TOLsplit, &
                                  MLMC_failprob, bnumMLMCsamps, num_ufunct, def_ufunct, spatial_Level, &
                                  num_benchsamps
  character(7) :: pltallopt                         !Plot all same opt
  character(3) :: default_ufunct    !default u functional or not?
  character(3) :: allL1s,allcenters !all basic L1/center based functionals?
  character(6) :: chosennorm        !norm to convert from text to number

  real(8)       :: dumreal,s2
  character(4)  :: setflags(3)
  character(20) :: dumchar !use this to "skip" a line

  integer :: i,default_ufunctc
  integer, allocatable :: tiarray2(:,:)

  open(unit=2,file="inputstoc.txt")

  read(2,*) rngseed
  read(2,*) probtype

  !--- Geometry - 'material' num of realz ---!
  read(2,*) dumchar
  read(2,*) numRealz,trannprt
  read(2,*) KLrnumRealz,KLrprintat

  !--- Geometry - 'material', 'Markov' type problem ---!
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

  !--- Geometry - 'material', Gauss or Gauss-based type problem ---!
  read(2,*) dumchar
  read(2,*) setflags(1),setflags(2)
  if(setflags(1)=='MB') flGBgeom = .false.
  if(setflags(2)=='same') flglGaussdiffrand = .false.
  read(2,*) setflags(1),chLNmode
  if(setflags(1)=='LN') flglLN = .true.
  read(2,*) GBsigave,GBsigvar
  read(2,*) GBscatrat
  read(2,*) GBlamc
  read(2,*) GBs

  !--- Large KL Options ---!
  read(2,*) dumchar
  read(2,*) KLres,KLrec
  read(2,*) numEigs
  read(2,*) binNumof

  !--- Large MCtrans Options ---!
  read(2,*) dumchar
  read(2,*) radMC,radWood,KLWood,LPMC,atmixMC,WAMC,GaussKL
  read(2,*) numParts
  read(2,*) LPamnumParts
  read(2,*) setflags(1)
  if(setflags(1)=='yes') flCorrMC    =.true.

  !--- Large MLMC Options ---!
  read(2,*) dumchar
  read(2,*) detMLMC
  read(2,*) MLMC_TOL, MLMC_TOLsplit

  !--- Lesser KL Options ---!
  read(2,*) dumchar
  read(2,*) setflags(1)
  if(setflags(1)=='tot') flmatbasedxs = .false.
  read(2,*) KLrnumpoints(2),KLrnumpoints(1)     !fixed xi, fixed point
  read(2,*) levsrefEig
  read(2,*) numrefinesameiter
  read(2,*) binSmallBound,binLargeBound
  read(2,*) KLnoise
  read(2,*) KLvarcalc,KLvarkept_tol
  read(2,*) numSlice
  read(2,*) Gaussrandtype

  !--- Lesser MCtrans Options ---!
  read(2,*) dumchar 
  read(2,*) rodOrplanar
  read(2,*) sourceType
  read(2,*) setflags(1),setflags(2)
  if(setflags(1)=='yes') flnegxs  =.true.
  if(setflags(2)=='yes') fldistneg=.true.
  read(2,*) setflags(1),meanadjust_tol
  if(setflags(1)=='yes') flmeanadjust=.true.
  read(2,*) refsigMode,userrefsig,negwgtbinnum,nwvalsperbin
  read(2,*) wgtmaxmin,wgtmax,wgtmin
  read(2,*) setflags(1),numLNxspts,numLNxsbins
  if(setflags(1)=='yes') flLNxscheck=.true.

  !--- Lesser MLMC Options ---!
  read(2,*) dumchar
  read(2,*) nextLevelFactor
  read(2,*) bnumMLMCsamps
  read(2,*) MLMC_failprob
  read(2,*) spatial_Level
  read(2,*) num_benchsamps

  !--- Geometry - 'coeffs' type problem ---!
  read(2,*) dumchar
  read(2,*) specialprob
  read(2,*) dumchar, nummat
  read(2,*) dumchar, s2
  if(probtype=='coeffs') s=s2
  read(2,*) dumchar, numcellsLevel0
  allocate(param1_mean(nummat))
  allocate(param1_uncert(nummat))
  allocate(param2_mean(nummat))
  allocate(param2_uncert(nummat))
  read(2,*) param1(1), param1_mean
  read(2,*) param1(2), param1_uncert
  read(2,*) param2(1), param2_mean
  read(2,*) param2(2), param2_uncert

  !--- Functionals - 'coeffs' type problem ---!
  read(2,*) dumchar
  read(2,*) num_ufunct
  read(2,*) default_ufunct,default_ufunctc
  read(2,*) allL1s,allcenters
  read(2,*) dumchar
  allocate(def_ufunct(num_ufunct,4))
  do i=1,num_ufunct  !read/set non-bulk values for functionals
    if(i==1 .and. default_ufunct=='yes') then
      def_ufunct(1,1) = 1               !first cell
      def_ufunct(1,2) = numcellsLevel0  !last cell
      def_ufunct(1,3) = 2               !L2 norm
      def_ufunct(1,4) = default_ufunctc !converge to this or not
    elseif(i==1) then
      read(2,*) def_ufunct(1,1),def_ufunct(1,2),chosennorm,def_ufunct(1,4)
      if(chosennorm=='L1')     def_ufunct(1,3)=1
      if(chosennorm=='L2')     def_ufunct(1,3)=2
      if(chosennorm=='center') def_ufunct(1,3)=3
    elseif(i>1) then
      read(2,*) def_ufunct(i,1),def_ufunct(i,2),chosennorm,def_ufunct(i,4)
      if(chosennorm=='L1')     def_ufunct(i,3)=1
      if(chosennorm=='L2')     def_ufunct(i,3)=2
      if(chosennorm=='center') def_ufunct(i,3)=3
    endif
  enddo

  if(allL1s=='yes') then  !add bulk values for functionals (many at one time)
    num_ufunct = num_ufunct+numcellsLevel0
    call move_alloc(def_ufunct,tiarray2)
    allocate(def_ufunct(num_ufunct,4))
    def_ufunct = 0
    def_ufunct(1:size(tiarray2(:,1)),:) = tiarray2
    do i=1,numcellsLevel0
      def_ufunct(size(tiarray2(:,1))+i,1) = i
      def_ufunct(size(tiarray2(:,1))+i,2) = i
      def_ufunct(size(tiarray2(:,1))+i,3) = 1
      def_ufunct(size(tiarray2(:,1))+i,4) = 0
    enddo
    deallocate(tiarray2)
  endif 
  if(allcenters=='yes') then
    num_ufunct = num_ufunct+numcellsLevel0
    call move_alloc(def_ufunct,tiarray2)
    allocate(def_ufunct(num_ufunct,4))
    def_ufunct = 0
    def_ufunct(1:size(tiarray2(:,1)),:) = tiarray2
    do i=1,numcellsLevel0
      def_ufunct(size(tiarray2(:,1))+i,1) = i
      def_ufunct(size(tiarray2(:,1))+i,2) = i
      def_ufunct(size(tiarray2(:,1))+i,3) = 3
      def_ufunct(size(tiarray2(:,1))+i,4) = 0
    enddo
    deallocate(tiarray2)
  endif 



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

  read(2,*) dumchar    !Plotting KLreconstructed realz
  read(2,*) pltKLrrealz(1),pltKLrrealz(2),pltKLrrealz(3),pltKLrrealz(4)
  read(2,*) pltKLrrealznumof
  allocate(pltKLrrealzwhich(2,pltKLrrealznumof))
  allocate(pltKLrrealzPointorXi(pltKLrrealznumof))
  do i=1,pltKLrrealznumof
    read(2,*) pltKLrrealzwhich(1,i),pltKLrrealzwhich(2,i),pltKLrrealzPointorXi(i)
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
  read(2,*) radMCbinplot,radWoodbinplot,KLWoodbinplot,GaussKLbinplot
  read(2,*) trprofile_binnum

  read(2,*) dumchar    !Plotting flux
  read(2,*) pltflux(1),pltflux(2),pltflux(3),pltflux(4)
  read(2,*) pltmatflux
  read(2,*) pltfluxtype
  read(2,*) fluxnumcells






  if( pltallopt .NE. 'default' ) then
    radMCbinplot   = pltallopt
    radWoodbinplot = pltallopt
    KLWoodbinplot  = pltallopt
    GaussKLbinplot = pltallopt
    pltEigf(1)     = pltallopt
    Corropts(1)    = pltallopt
    pltxiBins(1)   = pltallopt
    pltKLrrealz(1) = pltallopt
    pltgenrealz(1) = pltallopt
    pltCo(1)       = pltallopt
    pltflux(1)     = pltallopt
    pltmatflux     = pltallopt
  endif
  end subroutine readinputstoc










  subroutine testinputstoc
  use genRealzvars, only: sig, scatrat, numRealz, pltgenrealznumof, pltgenrealz, &
                          pltgenrealzwhich, flCorrMarkov, flCorrRealz
  use genSampvars, only: specialprob, param1, param2
  use KLvars, only: pltEigfwhich, pltxiBinswhich, pltCowhich, pltxiBinsnumof, pltEigfnumof, &
                    pltConumof, binNumof, numEigs, pltxiBins, pltEigf, pltCo, KLrnumpoints, &
                    KLrnumRealz, KLrprintat, pltKLrrealz, pltKLrrealznumof, pltKLrrealzwhich, &
                    pltKLrrealzPointorXi, KLres, KLrec, KLnoise, flmatbasedxs, &
                    Gaussrandtype, flCorrKL, flmeanadjust
  use MCvars, only: trannprt, sourceType, pltflux, radMC, radWood, KLWood, &
                    GaussKL, pltfluxtype, LPMC, atmixMC, radMCbinplot, radWoodbinplot, &
                    KLWoodbinplot, GaussKLbinplot, probtype, flnegxs
  use MLMCvars, only: def_ufunct, numcellsLevel0, nextLevelFactor
  integer :: fpointorxi(2)

  integer :: i, ones, zeros
  real(8) :: eps = 0.000001d0
  logical :: flstopstatus = .false., flsleep = .false.

  print *,"  "


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

                              !Test KLreconstruct print frequency
  if( (KLrprintat>KLrnumRealz .or. mod(KLrnumRealz,KLrprintat)/=0) .AND. KLrec=='yes' ) then 
    print *,"--Print to screen frequency for KLreconstruct must be factor of number of KLrealizations"
    flstopstatus = .true.
  endif

  fpointorxi = 0    !Test KLreconstruction num of realz, order of Eigs, and plot types
  do i=1,pltKLrrealznumof
    if( pltKLrrealzwhich(1,i)>KLrnumRealz .AND. pltKLrrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot more reconstructed realz than reconstructed"
      flstopstatus = .true.
    endif
    if( pltKLrrealzwhich(2,i)>numEigs .AND. pltKLrrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot reconstructed realz using more than calced num of Eigs"
      flstopstatus = .true.
    endif
    if( pltKLrrealzPointorXi(i) .NE. 'fpoint' .AND. pltKLrrealzPointorXi(i) .NE. 'fxi' ) then
      print *,"--User plot option for type of reconstruction needs to be either 'fpoint' or 'fxi'"
      flstopstatus = .true.
    endif
    !tally if more than one type present
    if(pltKLrrealzPointorXi(i) .EQ. 'fpoint') fpointorxi(1) = 1
    if(pltKLrrealzPointorXi(i) .EQ. 'fxi')    fpointorxi(2) = 1
  enddo
  !use tally to see if valid to plot
  if( fpointorxi(1) .NE. 0 .AND. fpointorxi(2) .NE. 0 .AND. pltKLrrealz(1) .NE. 'noplot' ) then
    if ( KLrnumpoints(1) .NE. KLrnumpoints(2) ) then
      print *,"--User must either plot only fpoint or fxi, or make num of points to plot same"
      flstopstatus = .true.
    endif
  endif
  if( KLres=='no' .and. KLrec=='yes' ) then
    KLres = 'yes'
    print *,"--User attempting KLrec w/o KLres, KLres has been set to 'yes'"
    flsleep = .true.
  endif
  if( GaussKL=='yes' .and. KLrec=='no' ) then
    KLrec = 'yes'
    print *,"--User attempting Gauss rand geom w/o KLrec, KLrec has been set to 'yes'"
    flsleep = .true.
  endif
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
  !Tests for Leakage pdf plotting options
  if(radMCbinplot  .ne. 'noplot'.or. radWoodbinplot .ne. 'noplot'.or. &
     KLWoodbinplot .ne. 'noplot'     ) then
    if(radMCbinplot .ne. 'noplot' .and. radMC=='no') then
      radMCbinplot = 'noplot'
      print *,"--User attempting to plot radMC leakage values w/o radMC, set to 'noplot'"
      flsleep = .true.
    endif
    if(radWoodbinplot .ne. 'noplot' .and. radWood=='no') then
      radWoodbinplot = 'noplot'
      print *,"--User attempting to plot radWood leakage values w/o radWood, set to 'noplot'"
      flsleep = .true.
    endif
    if(KLWoodbinplot .ne. 'noplot' .and. KLWood=='no') then
      KLWoodbinplot = 'noplot'
      print *,"--User attempting to plot KLWood leakage values w/o KLWood, set to 'noplot'"
      flsleep = .true.
    endif
    if(GaussKLbinplot .ne. 'noplot' .and. GaussKL=='no') then
      GaussKLbinplot = 'noplot'
      print *,"--User attempting to plot GaussKL leakage values w/o GaussKL, set to 'noplot'"
      flsleep = .true.
    endif
  endif

                              !Test radtransMC print frequency
  if( (trannprt>numRealz .or. mod(numRealz,trannprt)/=0) .AND. radMC=='yes' ) then 
    print *,"--Print to screen frequency for MCtran must be factor of number of realizations"
    flstopstatus = .true.
  endif

                              !Test plotting flux
  if( pltflux(1)/='noplot' .AND. radMC=='no' .AND. radWood=='no' .AND. KLWood=='no' &
      .and. LPMC=='no' .and. atmixMC=='no' .and. GaussKL=='no') then
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

  if( KLWood=='yes' ) then  !Tests for KLWood
    if( KLres=='no' .or. KLrec=='no' ) then
      !print *,"--User attempting to run KLWood w/o either KLresearch or KLreconstruct"
      !flstopstatus = .true.
      KLres='yes'
      KLrec='yes'
      print *,"--User attempting KLWood w/o either KLres or KLrec, both have been set to 'yes'"
      flsleep = .true.
    endif
    if( .not.flmatbasedxs .and. abs(scatrat(1)-scatrat(2))>eps ) then
      print *,"--User attempting to run KLWood w/ non-identical scattering ratios"
      flstopstatus = .true.
    endif
    if( KLWood=='yes' .and. numRealz/=KLrnumRealz ) then
      print *,"--User attempting to do transport over original and reconstructed with dif num of realz"
      flstopstatus = .true.
    endif
  endif
  
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

  if( KLnoise == 'yes' .AND. KLres == 'no' ) then !Test KLnoise w/o KLres
    print *,"--User trying to perform KLnoise without KLres"
    flstopstatus = .true.
  endif


  if( specialprob/='none' .and. specialprob/='mc2013.1' .and. specialprob/='mc2013.2' .and. &
      specialprob/='mc2015.1') then !Tests for genSampvars input
    print *,"--Enter 'none' or a valid special problem for genSamp specialproblem"
    flstopstatus = .true.
  endif
  if( param1(1)/='sigt' .and. param1(1)/='sigs' ) then
    print *,"--Enter 'sigt' or 'sigs' as first parameter type for genSamp"
    flstopstatus = .true.
  endif
  if( param2(1)/='c' .and. param2(1)/='siga' .and. param2(1)/='sigs' ) then
    print *,"--Enter 'c', 'siga', or 'sigs' as second parameter type for genSamp"
    flstopstatus = .true.
  endif
  if( (param1(1)=='sigt' .and. (param1(2)/='sigt1-abs' .and. param1(2)/='sigt1-frac')) .or. &
      (param1(1)=='sigs' .and. (param1(2)/='sigs1-abs' .and. param1(2)/='sigs1-frac')) ) then
    print *,"--Enter the uncertainty on param1 of a same type as the mean"
    flstopstatus = .true.
  endif
  if( (param2(1)=='c' .and. (param2(2)/='c1-abs' .and. param2(2)/='c1-frac')) .or. &
      (param2(1)=='siga' .and. (param2(2)/='siga1-abs' .and. param2(2)/='siga1-frac')) .or. &
      (param2(1)=='sigs' .and. (param2(2)/='sigs1-abs' .and. param2(2)/='sigs1-frac')) ) then
    print *,"--Enter the uncertainty on param2 of a same type as the mean"
    flstopstatus = .true.
  endif
  if( param1(1)=='sigs' .and. param2(1)/='siga' ) then
    print *,"--If param1 is sigs, param2 must be siga"
    flstopstatus = .true.
  endif 

  if(probtype=='coeffs') then
    ones  = 0
    zeros = 0
    do i=1,size(def_ufunct(:,1))
      if(def_ufunct(i,3)==3) then
        if(def_ufunct(i,1)/=def_ufunct(i,2)) then
          print *,"--Use 'center' option for functionals only over one cell"
          flstopstatus = .true.
        endif
        if(mod(nextLevelFactor,2)/=1) then
          print *,"--Use 'center' option only with odd refinement factor"
          flstopstatus = .true.
        endif
      endif
      if(def_ufunct(i,1)<1 .or. def_ufunct(i,2)>numcellsLevel0 .or. &
                                def_ufunct(i,1)>def_ufunct(i,2)) then
        print *,"--Define functionals over range of cells that exist"
        flstopstatus = .true.
      endif
      if(def_ufunct(i,4)==1) ones  = ones+1
      if(def_ufunct(i,4)==0) zeros = zeros+1
    enddo
    if(ones>1) then
      print *,"--Only first functional set to 1 will be required to converge"
      flsleep = .true.
    endif
    if(ones<1) then
      print *,"--Must choose a functional to converge"
      flstopstatus = .true.
    endif
    if(ones+zeros/=size(def_ufunct(:,1))) then
      print *,"--Set each functional's convergence status to either 0 or 1"
      flstopstatus = .true.
    endif
  endif

  if(flstopstatus) STOP 'killed'
  if(flsleep) call sleep(4)

  end subroutine testinputstoc



  subroutine global_allocate
  !This subroutine allocates and initializes all global variables
  use rngvars, only: rngappnum, rngseed
  use timevars, only: time, ntime, totparts, cumparts, FOM
  use genRealzvars, only: lam, P, s, numRealz, numPath, sumPath, sqrPath, largesti, &
                          totLength, lamc, sig, sigave, sigscatave, sigabsave, scatrat, &
                          flprint, numPosRealz, numNegRealz, numRealz, flGBgeom
  use KLvars, only: KLrrandarray, KLrnumpoints, numEigs, pltKLrrealznumof, KLrsig, &
                    KLrxisig, numSlice, gam, alpha, Ak, Eig, flMarkov, flGauss, &
                    xi, KLrxivals, KLrxivalss, pltKLrrealzarray, KLrnumRealz, flglGaussdiffrand, &
                    flGaussdiffrand
  use MCvars, only: fluxfaces, radMC, radWood, KLWood, WAMC, GaussKL, MCcaseson, MCcases, &
                    numParts, stocMC_reflection, stocMC_transmission, stocMC_absorption, &
                    numPosMCmeths, LPMC, atmixMC, LPamnumParts, stocMC_fluxall, &
                    stocMC_fluxmat1, stocMC_fluxmat2, pltflux, pltmatflux, &
                    fluxnumcells, flfluxplot
  use MLMCvars, only: MLMCcaseson, numPosMLMCmeths, MLMCcases, detMLMC
  use FEDiffSn, only: fliterstudy
  use mcnp_random, only: RN_init_problem
  integer :: i,icase

  !initialize rngvars
  rngappnum  = 0
  call RN_init_problem( 1, rngseed, int(0,8), int(0,8), 0)

  !allocate and initialize genRealzvars
  flprint    = .true.
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
  numPosRealz= 0
  numNegRealz= 0

  !allocate  KLresearch variables
  allocate(gam(numEigs))
  allocate(alpha(numEigs))
  allocate(Ak(numEigs))
  allocate(Eig(numEigs))
  allocate(xi(numRealz,numEigs))


  !allocate and initialize KLreconstruction variables
  allocate(KLrrandarray(KLrnumpoints(1),numEigs,pltKLrrealznumof+1))  !fpoint allocations
  allocate(KLrsig(KLrnumpoints(1)))
  allocate(KLrxivals(KLrnumRealz,numEigs))                            !fxi allocations
  if(GaussKL=='yes' .and. flglGaussdiffrand) allocate(KLrxivalss(KLrnumRealz,numEigs))
  flGaussdiffrand = .false.
  allocate(KLrxisig(KLrnumpoints(2)))
  allocate(pltKLrrealzarray(maxval(KLrnumpoints),pltKLrrealznumof+1)) !fpoint and/or fxi all


  !allocate/initialize MCvars
  allocate(MCcaseson(numPosMCmeths))
  MCcaseson = 0
  if(radMC  =='yes') MCcaseson(1) = 1
  if(radWood=='yes') MCcaseson(2) = 1
  if(KLWood =='yes') MCcaseson(3) = 1
  if(LPMC   =='yes') MCcaseson(4) = 1
  if(atmixMC=='yes') MCcaseson(5) = 1
  if(WAMC   =='yes') MCcaseson(6) = 1
  if(GaussKL=='yes') MCcaseson(7) = 1


  allocate(MCcases(numPosMCmeths))
  MCcases(1) = 'radMC'
  MCcases(2) = 'radWood'
  MCcases(3) = 'KLWood'
  MCcases(4) = 'LPMC'
  MCcases(5) = 'atmixMC'
  MCcases(6) = 'WAMC'
  MCcases(7) = 'GaussKL'

  allocate(stocMC_reflection(numPosMCmeths,2))   !global MC variables for each method
  allocate(stocMC_transmission(numPosMCmeths,2)) !rank 2 holds 1=average, 2=deviation
  allocate(stocMC_absorption(numPosMCmeths,2))
  stocMC_reflection   = 0.0d0
  stocMC_transmission = 0.0d0
  stocMC_absorption   = 0.0d0

  if(sum(MCcaseson(1:7))>0 .or. (.not.flGBgeom .and. MCcaseson(7)==1)) flMarkov = .true.
  if(flGBgeom .and. MCcaseson(7)==1) flGauss = .true.
  
  flfluxplot = .false.  !flux variable allocations
  if( pltflux(1)=='plot' .or. pltflux(1)=='preview' .or. &
      pltmatflux=='plot' .or. pltmatflux=='preview' ) flfluxplot = .true.
  if(flfluxplot) then
    allocate(fluxfaces(fluxnumcells+1))
    fluxfaces = 0.0d0
    do i=1,fluxnumcells+1
      fluxfaces(i) = (s/fluxnumcells) * (i-1)
    enddo
  endif
  if( pltflux(1)=='plot' .or. pltflux(1)=='preview' .or. &!mat irrespective flux allocations
    (flfluxplot .and. (KLWood=='yes' .or. atmixMC=='yes')) ) then !KLWood, atmixMC, respective stored
                                                                  !here (actually irresective)
    allocate(stocMC_fluxall(fluxnumcells,numPosMCmeths,2))
    stocMC_fluxall = 0.0d0
  endif
  if( pltmatflux=='plot' .or. pltmatflux=='preview' ) then !mat respective flux allocations
    allocate(stocMC_fluxmat1(fluxnumcells,numPosMCmeths,2))
    allocate(stocMC_fluxmat2(fluxnumcells,numPosMCmeths,2))
    stocMC_fluxmat1 = 0.0d0
    stocMC_fluxmat2 = 0.0d0
  endif


  !allocate/initialize MLMCvars
  allocate(MLMCcaseson(numPosMLMCmeths))
  MLMCcaseson = 0
  if(detMLMC=='detMLMC') MLMCcaseson(1) = 1
  if(detMLMC=='spatial') MLMCcaseson(2) = 1
  if(detMLMC=='iter')    MLMCcaseson(3) = 1
  if(detMLMC=='bench')   MLMCcaseson(4) = 1

  allocate(MLMCcases(numPosMLMCmeths))
  MLMCcases(1) = 'detMLMC'
  MLMCcases(2) = 'spatial'
  MLMCcases(3) = 'iter'
  MLMCcases(4) = 'bench'

  !initialize FEDiffSn variables
  fliterstudy = .false.

  !allocate and initialize timevars
  allocate(time(ntime))
  time = 0.0d0
  allocate(FOM(ntime,2))
  FOM  = 0.0d0
  allocate(totparts(numPosMCmeths))
  allocate(cumparts(numPosMCmeths))
  totparts = 0
  cumparts = 0
  do icase=1,numPosMCmeths
    if(MCcaseson(icase)==1) then
      if(MCcases(icase)=='LPMC' .or. MCcases(icase)=='atmixMC') then
        totparts(icase) = LPamnumParts
      else
        totparts(icase) = numRealz*numParts
      endif
    endif
  enddo

  end subroutine global_allocate



  subroutine clearreports
  call system("rm texts/*.out")
  end subroutine clearreports



  subroutine finalreport
  call system("cat texts/Acase.out texts/Woodnegstats.out texts/MCleakage.out texts/timereport.out > texts/finalreport.out")
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








end module Loadcase

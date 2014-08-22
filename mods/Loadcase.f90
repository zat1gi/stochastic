module Loadcase
  implicit none

CONTAINS
  ! print statements in this module use # 100-199

  subroutine Acase_load
  use genRealzvars, only:   Adamscase, sig, lam, scatrat, s  

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

  end subroutine Acase_load





  subroutine Acase_print
  use genRealzvars, only: Adamscase, sig, lam, scatrat, s, lamc

    lamc = lam(1)*lam(2)/(lam(1)+lam(2))

    100 format("  Values loaded for special Adamscase :",f5.1)
    101 format("    sig(1)    :",f10.6,"   sig(2)    :",f10.6)
    102 format("    scatrat(1):",f6.2,"       scatrat(2):",f6.2)
    103 format("    lam(1)    :",f6.2,"       lam(2)    :",f6.2)
    104 format("    s         :",f6.2,"       lamc      :",f6.2)

    write(*,*)
    write(*,100) Adamscase
    write(*,101) sig(1),sig(2)
    write(*,102) scatrat(1),scatrat(2)
    write(*,103) lam(1),lam(2)
    write(*,104) s,lamc
    write(*,*)

  end subroutine Acase_print





  subroutine readinputstoc( KLres,KLnoise,&
                            KLrec,&
                            numParts,trannprt,radMC,rodOrplanar,results,&
                            radWood,KLWood,allowneg,&
                            distneg,plotflux,pfnumcells,pltflux,sourceType,seed )
  use genRealzvars,         only: Adamscase, sig, scatrat, lam, s, numRealz, pltgenrealznumof, &
                                  pltgenrealz, plotmatdxs, pltgenrealzwhich
  use KLvars,               only: KLvarcalc, KLvarkept_tol, pltEigfwhich, pltxiBinswhich, &
                                  pltCowhich, pltxiBinsnumof, pltEigfnumof, pltConumof, binNumof,&
                                  numEigs, numSlice, levsrefEig, Corrnumpoints, binSmallBound, &
                                  binLargeBound, pltxiBins, pltxiBinsgauss, pltEigf, pltCo, &
                                  Corropts, KLrnumpoints, KLrnumRealz, KLrprintat, pltKLrrealz, &
                                  pltKLrrealznumof, pltKLrrealzwhich, pltKLrrealzPointorXi
  use MCvars,               only: trprofile_binnum, radMCbinplot, radWoodbinplot, KLWoodbinplot
  use KLmeanadjust,         only: KLadjust, meanadjust_tol
  integer :: seed                                   !adv seed
  character(3) :: KLres,KLnoise
  character(3) :: KLrec
  integer :: numParts,trannprt                      !radtransMC opts
  character(3) :: radMC
  character(6) :: rodOrplanar,results,sourceType
  character(3) :: radWood,KLWood,allowneg,distneg   !Woodcock opts
  integer :: pfnumcells                             !rad & Wood flux plotting
  character(6) :: plotflux(2)
  character(7) :: pltflux(4)

  character(7) :: pltallopt                         !Plot all same opt

  real(8)       :: dumreal
  character(20) :: dumchar !use this to "skip" a line

  integer :: i

  open(unit=2,file="inputstoc.txt")

  !--- genRealz ---!
  read(2,*) dumchar    !-- genRealz Shared Variables --!
  read(2,*) Adamscase
  read(2,*) sig(1),sig(2)
  read(2,*) scatrat(1),scatrat(2)
  read(2,*) lam(1),lam(2)
  read(2,*) s
  read(2,*) numRealz


  !--- KL ---!
  read(2,*) dumchar    !-- KL Shared Variables --!
  read(2,*) KLvarcalc
  read(2,*) KLvarkept_tol


  !--- Transport ---!
  read(2,*) dumchar    !-- Transport Shared Variables --!
  read(2,*) trprofile_binnum

  read(2,*) dumchar    !-- radMC (TMC) Specific Variables --!
  read(2,*) dumchar      !!-radMCbinplot
  read(2,*) radMCbinplot

  read(2,*) dumchar    !-- radWood (WMC) Specific Variables --!
  read(2,*) dumchar      !!-radWoodbinplot
  read(2,*) radWoodbinplot

  read(2,*) dumchar    !-- KLWood (KLWMC) Specific Variables --!
  read(2,*) dumchar      !!-KLWoodbinplot
  read(2,*) KLWoodbinplot




  read(2,*) seed

  read(2,*) dumchar    !KL research options
  read(2,*) KLres
  read(2,*) binNumof
  read(2,*) numEigs
  read(2,*) numSlice
  read(2,*) levsrefEig
  read(2,*) binSmallBound,binLargeBound
  read(2,*) KLnoise

  read(2,*) dumchar    !KL reconstruct options
  read(2,*) KLrec
  read(2,*) KLrnumpoints(1)     !fixed point
  read(2,*) KLrnumpoints(2)     !fixed xi
  read(2,*) KLrnumRealz
  read(2,*) KLrprintat

  read(2,*) dumchar    !radtransMC options
  read(2,*) radMC
  read(2,*) numParts
  read(2,*) trannprt
  read(2,*) rodOrplanar
  read(2,*) sourceType
  read(2,*) results

  read(2,*) dumchar    !Woodcock options
  read(2,*) radWood
  read(2,*) KLWood
  read(2,*) allowneg,distneg
  read(2,*) KLadjust,meanadjust_tol


  read(2,*) dumchar    !All Plot Same Way Option
  read(2,*) pltallopt

  read(2,*) dumchar    !Plotting flux
  read(2,*) pltflux(1),pltflux(2),pltflux(3),pltflux(4)
  read(2,*) plotmatdxs
  read(2,*) plotflux(1),plotflux(2)
  read(2,*) pfnumcells

  read(2,*) dumchar    !Plotting Eigenfunction
  read(2,*) pltEigf(1),pltEigf(2),pltEigf(3),pltEigf(4)
  read(2,*) pltEigfnumof
  allocate(pltEigfwhich(pltEigfnumof))
  do i=1,pltEigfnumof
    read(2,*) pltEigfwhich(i)
  enddo

  read(2,*) dumchar    !Plotting Correlation contours
  read(2,*) Corropts(1),Corropts(2)
  read(2,*) Corrnumpoints

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

  read(2,*) dumchar    !Plotting genRealz realz
  read(2,*) pltgenrealz(1),pltgenrealz(2),pltgenrealz(3),pltgenrealz(4)
  read(2,*) pltgenrealznumof
  allocate(pltgenrealzwhich(pltgenrealznumof))
  do i=1,pltgenrealznumof
    read(2,*) pltgenrealzwhich(i)
  enddo

  read(2,*) dumchar    !Plotting Variace (Co)
  read(2,*) pltCo(1),pltCo(2),pltCo(3),pltCo(4)
  read(2,*) pltConumof
  allocate(pltCowhich(2,pltConumof))
  do i=1,pltConumof
    read(2,*) pltCowhich(1,i),pltCowhich(2,i)
  enddo




  if( pltallopt .NE. 'default' ) then
    pltEigf(1)     =pltallopt
    Corropts(1)    =pltallopt
    pltxiBins(1)   =pltallopt
    pltKLrrealz(1) =pltallopt
    pltgenrealz(1) =pltallopt
    pltCo(1)       =pltallopt
    pltflux(1)     =pltallopt
  endif

  end subroutine readinputstoc










  subroutine testinputstoc( trannprt,KLres,KLrec,radWood,&
                            radMC,&
                            KLnoise,KLWood,pltflux,&
                            sourceType,allowneg,distneg )
  use genRealzvars, only: sig, scatrat, numRealz, pltgenrealznumof, pltgenrealz, &
                          pltgenrealzwhich
  use KLvars, only: pltEigfwhich, pltxiBinswhich, pltCowhich, pltxiBinsnumof, pltEigfnumof, &
                    pltConumof, binNumof, numEigs, pltxiBins, pltEigf, pltCo, KLrnumpoints, &
                    KLrnumRealz, KLrprintat, pltKLrrealz, pltKLrrealznumof, pltKLrrealzwhich, &
                    pltKLrrealzPointorXi
  integer :: trannprt
  integer :: fpointorxi(2)
  character(6) :: sourceType
  character(7) :: pltflux(4)
  character(3) :: KLres,KLrec,radMC,KLnoise,KLWood,radWood,allowneg,distneg

  integer :: i
  real(8) :: smallersig,largersig,sigratio
  character(3) :: stopstatus = 'no', run = 'no'

  print *,"  "


  do i=1,pltEigfnumof    !Test Eigenfunction plotting order of Eigs
    if( pltEigfwhich(i)>numEigs .AND. pltEigf(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot Eigenfunction of higher order than Eigenvalues calculated"
      stopstatus = 'yes'
    endif
  enddo

  do i=1,pltxiBinsnumof  !Test xiBins plotting order of Eigs and num of bins
    if( pltxiBinswhich(1,i)>numEigs .AND. pltxiBins(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot xiBins for Eigenvalue of higher order than calculated"
      stopstatus = 'yes'
    endif
    if( pltxiBinswhich(2,i)>numRealz .AND. pltxiBins(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot xibins for more realizations than generated"
      stopstatus = 'yes'
    endif
  enddo

  if( KLnoise == 'yes' .AND. KLres == 'no' ) then !Test KLnoise w/o KLres
    print *,"--User tryint to perform KLnoise without KLres"
    stopstatus = 'yes'
  endif

  if( KLrprintat>KLrnumRealz ) then  !Test print frequency on reconstruction
    print *,"--User trying to print KLreconstruction update less often than total reconstructions"
    stopstatus = 'yes'
  endif

  fpointorxi = 0    !Test KLreconstruction num of realz, order of Eigs, and plot types
  do i=1,pltKLrrealznumof
    if( pltKLrrealzwhich(1,i)>KLrnumRealz .AND. pltKLrrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot more reconstructed realz than reconstructed"
      stopstatus = 'yes'
    endif
    if( pltKLrrealzwhich(2,i)>numEigs .AND. pltKLrrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot reconstructed realz using more than calced num of Eigs"
      stopstatus = 'yes'
    endif
    if( pltKLrrealzPointorXi(i) .NE. 'fpoint' .AND. pltKLrrealzPointorXi(i) .NE. 'fxi' ) then
      print *,"--User plot option for type of reconstruction needs to be either 'fpoint' or 'fxi'"
      stopstatus = 'yes'
    endif
    !tally if more than one type present
    if(pltKLrrealzPointorXi(i) .EQ. 'fpoint') fpointorxi(1) = 1
    if(pltKLrrealzPointorXi(i) .EQ. 'fxi')    fpointorxi(2) = 1
  enddo
  !use tally to see if valid to plot
  if( fpointorxi(1) .NE. 0 .AND. fpointorxi(2) .NE. 0 .AND. pltKLrrealz(1) .NE. 'noplot' ) then
    if ( KLrnumpoints(1) .NE. KLrnumpoints(2) ) then
      print *,"--User must either plot only fpoint or fxi, or make num of points to plot same"
      stopstatus = 'yes'
    endif
  endif
  if( KLres=='no' .AND. KLrec=='yes' ) then
    print *,"--User attempting to reconstuct realizations (KLrec) without harvesting values (KLres)"
    print *,"--KLrec has been switched to 'no'"
    KLrec = 'no'
  endif

  if( trannprt>numRealz .AND. radMC=='yes' ) then !Test radtransMC print frequency
    print *,"--User attempting to print to screen after more realizations than calculated (radtrans)"
    stopstatus = 'yes'
  endif

                              !Test plotting flux
  if( pltflux(1)/='noplot' .AND. radMC=='no' .AND. radWood=='no' .AND. KLWood=='no' ) then
    print *,"--User attempting to plot flux when no transport calculations are made"
    stopstatus = 'yes'
  endif
  if( sourceType/='left' .AND. sourceType/='intern' ) then
    print *,"--User attempting to run invalid source type.  Please put either 'left' or 'intern'"
    stopstatus = 'yes'
  endif

  do i=1,pltgenrealznumof    !Test genRealz plotting over selected realz
    if( pltgenrealzwhich(i)>numRealz .AND. pltgenrealz(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot realizations that are not created"
      stopstatus = 'yes'
    endif
  enddo

  do i=1,pltConumof          !Test plotCo for Eig choice and CoEffExp or CoEffAct option
    if( pltCowhich(1,i)>numEigs .AND. pltCo(1) .NE. 'noplot' ) then
      print *,"--User attempting to plot CoEff values for eigenvalues not calculated"
      stopstatus = 'yes'
    endif
    if( pltCowhich(2,i)/=1 .AND. pltCowhich(2,i)/=2 .AND. pltCo(1) .NE. 'noplot' ) then
      print *,"--User input for 'CoEffExp vs CoEffAct' not valid.  Enter a '1' or a '2'"
      stopstatus = 'yes'
    endif
  enddo

  if( KLWood=='yes' ) then  !Tests for KLWood
    if( KLres=='no' .OR. KLrec=='no' ) then
      print *,"--User attempting to run KLWood w/o either KLresearch or KLreconstruct"
      stopstatus = 'yes'
    endif
    if( scatrat(1) /= scatrat(2) ) then
      print *,"--User attempting to run KLWood w/ non-identical scattering ratios"
      stopstatus = 'yes'
    endif
    smallersig = minval(sig)
    largersig  = maxval(sig)
    sigratio   = (largersig-smallersig)/smallersig
    if( sigratio > 0.33334d0 .AND. stopstatus=='no' .and. allowneg=='no') then
      print *,"--User attempting to run KLWood where neg reconstructed xs values may exist"
      print *,"   -if you choose to run this, you will want your # of pnts to recon at to be quite high"
      print *,"   -please either 'run' to run anyway, or anything else to exit"
      read(*,*) run
      if( run .NE. 'run' ) stopstatus = 'yes'
    endif
    if( allowneg=='no' .and. distneg=='yes' ) then
      print *,"--User attempting to redistribute negative xs values without allowneg on"
      stopstatus = 'yes'
    endif
    if( radMC=='yes' .and. numRealz/=KLrnumRealz ) then
      print *,"--User attempting to do transport over original and reconstructed with dif num of realz"
      stopstatus = 'yes'
    endif
  endif


  if( stopstatus=='yes' ) STOP 'killed'

  end subroutine testinputstoc








  subroutine timereport( runtime,time,ntime,KLres,KLrec,radMC,radWood,KLWood,&
                         KLnoise )
  integer :: ntime
  real(8) :: runtime,time(:)
  character(3) :: KLres,KLrec,radMC,radWood,KLWood,KLnoise

  integer :: i
  real(8) :: othertime,otherpercent
  real(8),allocatable :: pertime(:)

  allocate(pertime(ntime))


  othertime = runtime
  do i=1,ntime
    time(i)    = time(i) / 60
    othertime  = othertime - time(i)
    pertime(i) = time(i) / runtime * 100
  enddo
  otherpercent = othertime / runtime * 100

  110 format("    genRealz   :",f10.2," min     per:",f6.2," %")
  111 format("    radMC      :",f10.2," min     per:",f6.2," %")
  112 format("    radWood    :",f10.2," min     per:",f6.2," %")
  113 format("    KLnoise    :",f10.2," min     per:",f6.2," %")
  114 format("    KLcol      :",f10.2," min     per:",f6.2," %")
  115 format("    KLrec      :",f10.2," min     per:",f6.2," %")
  116 format("    KLWood     :",f10.2," min     per:",f6.2," %")
  117 format("    other      :",f10.2," min     per:",f6.2," %")
!  118 format("    total      :",f10.2," min")

                     write(*,110) time(1),pertime(1)
  if(radMC=='yes')   write(*,111) time(2),pertime(2)
  if(radWood=='yes') write(*,112) time(3),pertime(3)
  if(KLnoise=='yes') write(*,113) time(4),pertime(4)
  if(KLres=='yes')   write(*,114) time(5),pertime(5)
  if(KLrec=='yes')   write(*,115) time(6),pertime(6)
  if(KLWood=='yes')  write(*,116) time(7),pertime(7)
                     write(*,117) othertime,otherpercent
!  write(*,118) runtime


  end subroutine timereport









end module Loadcase

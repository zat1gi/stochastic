module KLresearch
  use utilities
  use KLmeanadjust
  use timeman
  implicit none

CONTAINS
  ! print statements in this module use # 400-499






  subroutine KL_eigenvalue( numEigs,P,sigave,&
                            levsrefEig,lamc,numSlice,pltEigf,&
                            KLrxivals,KLrnumRealz )
  !This subroutine: 1) calcaltes some initial values used here and later
  !2) Solves the transcendental equation which yields gamma
  !3) From gamma solves: alpha, lambda (Eigenvalue), & the normalization const A_k
  !4) Prints and plots Eigenfunctions if input specifies
  !5) Calculates variance maintained with # of eigvals if input specifies
  use genRealzvars, only: sig, lam, s, numRealz
  use KLvars,       only: KLvarkept_tol, KLvarcalc, AllEig, Allgam, varmain, gam, alpha,&
                          Ak, Eig, xi, pltEigfwhich, pltEigfnumof
  integer :: numEigs,levsrefEig,numSlice,KLrnumRealz
  real(8),allocatable :: KLrxivals(:,:)
  real(8) :: P(2),sigave,lamc
  character(7) :: pltEigf(4)

  real(8) :: stepGam=0 !if 0 code chooses
  integer :: index,l,level,curEig,i,j
  real(8) :: refstepGam,TT,curGam,Co,sqrtEig(numEigs)
  real(8) :: absdiff,absdiff_1=0,absdiff_2=0,testval(numEigs),sliceSize
  real(8),allocatable :: Eigfplotarray(:,:)
  integer :: prevsize,newsize
  real(8) :: Eigval,Eigvalsum,oldsEig,newsEig
  real(8), allocatable :: oldAllgam(:),oldAllEig(:),cumeig(:)

  allocate(Eigfplotarray(numSlice,numSlice+1))
  allocate(gam(numEigs))
  allocate(alpha(numEigs))
  allocate(Ak(numEigs))
  allocate(Eig(numEigs))
  allocate(xi(numRealz,numEigs))
  allocate(KLrxivals(KLrnumRealz,numEigs))

  424 format("   P(1):",f8.5,"   P(2)  :",f8.5,"   lamc:",f8.5,"   TT:",f8.5)
  425 format("   Co  :",f8.5,"   sigave:",f8.5,"   stepGam:",f8.5)
  lamc   =(lam(1)*lam(2))/(lam(1)+lam(2))
  P(1)     =lam(1)/(lam(1)+lam(2))
  P(2)     =lam(2)/(lam(1)+lam(2))
  TT       =s/lamc
  sigave   =P(1)*sig(1)+P(2)*sig(2)
  Co       =P(1)*P(2)*(sig(1)-sig(2))**2
  if( stepGam==0 ) then
    stepGam  =1/TT/50
  endif

  !Initial guesses for Gam
  420 format(i7,"     ",f15.9," ",f19.14)

  allocate(AllEig(10))
  allocate(Allgam(10))
  Allgam    = 0d0
  prevsize  = 1
  newsize   = 10
  Eigvalsum = 0d0
  write(*,*) " Calculating Eigenvalues for computation and varaince "
  write(*,*) "          Eigindx    Gam vals       Eig vals         tol check         tol "
  460 format("Eigcalc: "i7,"   ",f12.6,"   ",f12.9,"     ",f12.9,"    ",f12.9)

  do     !until tolerance
    !write(*,*)
    !write(*,*) " Finding apprx Gam vals for Eigenvalue calcs"
    !write(*,*) "  Eigindx    apprx Gam vals        absdiff"
    curGam=Allgam(prevsize)+stepGam
    curEig=prevsize
    do while ( curEig<=newsize )

      absdiff=abs( tan(curGam*TT)-(2d0*curGam)/(curGam**2-1d0) )
      open(unit=11, file="absdiffGam.txt")
      write(11,420) curEig,curGam,absdiff

      if( absdiff_2>absdiff_1 .AND. absdiff>absdiff_1 ) then
        !write(*,420) curEig,curGam,absdiff_1
        Allgam(curEig)=curGam
        curEig=curEig+1
      endif

      absdiff_2=absdiff_1
      absdiff_1=absdiff
      curGam=curGam+stepGam
    enddo
    close(unit=11)

    !Refine guesses for gam within tolerance
    !write(*,*)
    !write(*,*) " Refining Gamma values for Eigenvalue calcs"
    !write(*,*) "  Eigindx    apprx Gam vals        absdiff          refine lev     inttest=?1"

    !421 format(i7,"     ",f15.9,"     ",f15.14,"      Level:",i3)

    do curEig=prevsize,newsize                  !loop each Eigenvalue
      refstepGam=stepGam

      do level=1,levsrefEig                 !loop levels of refinement
        curGam=Allgam(curEig)-11d0*refstepGam
        refstepGam=refstepGam/10
        absdiff_1=0
        absdiff_2=0

        do l=1,220                            !loop through range & test
          absdiff=abs( tan(curGam*TT)-(2d0*curGam)/(curGam**2-1d0) )
          if( absdiff_2>absdiff_1 .AND. absdiff>absdiff_1 ) then
            !write(*,421) curEig,curGam,absdiff_1,level
            !if(level==levsrefEig) write(*,421) curEig,curGam,absdiff_1,level
            Allgam(curEig)=curGam
          endif
          absdiff_2=absdiff_1
          absdiff_1=absdiff
          curGam=curGam+refstepGam
        enddo
      enddo
      AllEig(curEig)     =2d0*Co*lamc/(Allgam(curEig)**2+1d0)
      Eigvalsum = Eigvalsum + AllEig(curEig)
    enddo
    Eigval = AllEig(newsize)

    !end loop?
    oldsEig = sqrt(Eigvalsum-Eigval)
    newsEig = sqrt(Eigvalsum)
    write(*,460) newsize,Allgam(newsize),AllEig(newsize),abs(oldsEig-newsEig)/newsEig,KLvarkept_tol
    if(KLvarcalc=='yes' .and. abs(oldsEig-newsEig)/newsEig<KLvarkept_tol .and. newsize>numEigs) exit
    if(KLvarcalc=='no'  .and. newsize>=numEigs) exit

    !tack 10 new array values on
    if(.not. allocated(oldAllEig)) allocate(oldAllEig(newsize))
    if(.not. allocated(oldAllgam)) allocate(oldAllgam(newsize))
    oldAllEig = AllEig
    oldAllgam = Allgam
    deallocate(AllEig)
    deallocate(Allgam)
    prevsize = prevsize+10
    newsize  = newsize +10
    allocate(AllEig(newsize))
    allocate(Allgam(newsize))
    AllEig(1:prevsize-1) = oldAllEig(1:prevsize-1)
    Allgam(1:prevsize-1) = oldAllgam(1:prevsize-1)
    Allgam(prevsize) = oldAllgam(prevsize-1) !for curGam starting position in next round
    deallocate(oldAllEig)
    deallocate(oldAllgam)    
  enddo
  call system("mv absdiffGam.txt texts")

  Eig(1:numEigs) = AllEig(1:numEigs)
  gam(1:numEigs) = Allgam(1:numEigs)

  write(*,*) "     Int test"
  do curEig=1,numEigs
    !Calc other values like alpha, norm const (Ak), eigenvalue, etc.
    alpha(curEig)   =gam(curEig)/lamc
    Ak(curEig)      =sqrt(1d0/(  s/2d0*(gam(curEig)**2+1d0)+lamc  ))
    !Eig(curEig)     =2d0*Co*lamc/(gam(curEig)**2+1d0)
    sqrtEig(curEig) =sqrt(Eig(curEig))
    !integrate to 1 tests
    427 format("  ",f13.7,"   Ak:",f13.7)
    testval(curEig)=Ak(curEig)**2*0.5d0*( s*(1+gam(curEig)**2)+&
                    lamc*(2d0- (sin(alpha(curEig)*s))**2 ) )
    write(*,427) testval(curEig),Ak(curEig)

    !422 format("  curEig:",i4,"   gam:",f16.7,"   alpha:",f16.7)
    !423 format("                 Ak:",f16.7,"   Eig  :",f16.7)
    !write(*,422) curEig,gam(curEig),alpha(curEig)
    !write(*,423) Ak(curEig),Eig(curEig)
    !write(*,*)
  enddo

  allocate(cumeig(newsize))
  allocate(varmain(newsize))
  do curEig=1,newsize
    cumeig(curEig) = sum(AllEig(1:curEig))
  enddo
  varmain = sqrt(cumeig)/sqrt(cumeig(newsize))


  write(*,*)
  if(KLvarcalc=='no')  write(*,*) "    Eigenvalues and their contributions"
  if(KLvarcalc=='no')  write(*,*) "  Eigindx       Eig vals        sqrt(Eig)"
  if(KLvarcalc=='yes') write(*,*) "    Eigenvalues, their contributions, and KL maintained variance"
  if(KLvarcalc=='yes') write(*,*) "  Eigindx       Eig vals        sqrt(Eig)       var maint"
  426 format(i7,"     ",f15.9,"  ",f13.9,"     ",f13.9,"   ",f13.9)
  do curEig=1,numEigs
    if(KLvarcalc=='no')  write(*,426) curEig,Eig(curEig),sqrtEig(curEig)
    if(KLvarcalc=='yes') write(*,426) curEig,Eig(curEig),sqrtEig(curEig),varmain(curEig)
  enddo

  428 format("   lamc:           ",f8.3)
  429 format("   TT  :           ",f13.8)
  write(*,428) lamc
  write(*,429) TT
  write(*,*)

  !This section plots the Eigenfunctions according to the options in the input file
  if( pltEigf(1) .NE. 'noplot' ) then  !plot using generic plotter
    sliceSize = s/numSlice
    do i=1,numSlice
      Eigfplotarray(i,1) = sliceSize*i-sliceSize/2
      do j=1,pltEigfnumof
        curEig=pltEigfwhich(j)
        Eigfplotarray(i,j+1) = Eigfunc(Ak(curEig),alpha(curEig),lamc,Eigfplotarray(i,1))

      enddo
    enddo

    call generic_plotter( numSlice,pltEigfnumof,Eigfplotarray,pltEigf )

    call system("mv genericplot.txt plots/Eigfplot/Eigfplot.txt")
    call system("mv genericplot.ps  plots/Eigfplot/Eigfplot.ps")
    call system("mv genericplot.pdf plots/Eigfplot/Eigfplot.pdf")
  endif


  xi = 0d0  !initializing xi=0 for all values

  end subroutine KL_eigenvalue







  subroutine KL_Correlation( Corropts,Corrnumpoints,numEigs,&
                             lamc,sigave,CoExp,P )
  !This subroutine calculates Correlation between two points in a
  !realization based upon the expected value, and the observed 
  !value (function of Eigenfunctions and values).
  !It then plots in 3D if user has specified.
  use genRealzvars, only: sig, s
  use KLvars, only: alpha, Ak, Eig
  integer :: Corrnumpoints,numEigs
  real(8) :: lamc,sigave,CoExp,P(2)
  character(7) :: Corropts(2)

  integer :: x,y,curEig
  real(8) :: stepsize,curx,cury,Eigfx,Eigfy
  real(8) :: Correxpect(Corrnumpoints,Corrnumpoints)
  real(8) :: Corryield(Corrnumpoints,Corrnumpoints)
  character(17) :: Corrchar = '"Correlation.txt"' 
  character(22) :: Corre    = '"Correlation Expected"'
  character(22) :: Corry    = '"Correlation Yielded "'

  CoExp= P(1)*P(2) * (sig(1)-sig(2))**2 !first calc of vairance

  if( Corropts(1) .NE. "noplot" ) then  !only perform if "plot" or "preview"

  Correxpect = 0  !initialize
  Corryield  = 0

  stepsize = s/(Corrnumpoints)  !set up stepsize


  do x=1,Corrnumpoints  !cycle through x and y
    if(x==1) curx=stepsize/2 
    do y=1,Corrnumpoints
      if(y==1) cury=stepsize/2

        !calc Correxpect
      Correxpect(x,y) = CoExp * exp( - abs(curx-cury)/lamc )

        !sum Eigs to calc Corryield
      do curEig=1,numEigs
        Eigfx = Eigfunc(Ak(curEig),alpha(curEig),lamc,curx)
        Eigfy = Eigfunc(Ak(curEig),alpha(curEig),lamc,cury)
        Corryield(x,y) = Corryield(x,y) + Eig(curEig) * Eigfx * Eigfy
      enddo

        !print Correxpect and Corryield
      210 format("#       x                y            Correxpect        Corryield")
      211 format(" ",f12.6,"     ",f12.6,"     ",f12.6,"     ",f12.6)
      open(unit=15,file="Correlation.txt")
      if(x==1 .AND. y==1) write(15,210)
      write(15,211) curx,cury,Correxpect(x,y),Corryield(x,y)

      cury=cury+stepsize
    enddo !y-cycle
    curx=curx+stepsize
  enddo !x-cycle
  close(unit=15) !close file


  open(unit=16,file="Correlationfirst.gnu") !make custom gnu script
  212 format("set dgrid3d ",i10,",",i10)
  213 format("splot ",A17," using 1:2:3 with lines title ",A22)
  214 format("splot ",A17," using 1:2:4 with lines title ",A22)
  215 format("splot ",A17," using 1:2:3 with lines title ",A22,", ",A17," using 1:2:4 with lines title ",A22)
  219 format(A40)
  write(16,212) Corrnumpoints,Corrnumpoints
  write(16,219)'set view 60,130                         '
  write(16,219)'                                        '
  write(16,219)'set title "Positional Correlation Plots"'
  write(16,219)'set ylabel "\"y\" or \"x2\" position"   '
  write(16,219)'set xlabel "\"x\" or \"x1\" position"   '
  write(16,219)'                                        '
  if(Corropts(2)=='expect')  write(16,213) Corrchar,Corre
  if(Corropts(2)=='yield')   write(16,214) Corrchar,Corry
  if(Corropts(2)=='both')    write(16,215) Corrchar,Corre,Corrchar,Corry
  close(unit=16)

  if(Corropts(1)=='preview') then
    call system("cat Correlationfirst.gnu plots/Correlationgnus/Corrps_pause.txt > Correlation.gnu")
  else
    call system("cat Correlationfirst.gnu plots/Correlationgnus/Corrps.txt > Correlation.gnu")
  endif

  call system("gnuplot Correlation.gnu") !plot using gnuplot

  call system("ps2pdf Correlation.ps Correlation.pdf")         !convert to pdf
  call system("ps2eps Correlation.ps")         !convert to eps

  call system("rm Correlationfirst.gnu")                       !clean up space
  call system("rm plots/Correlation.eps")
  call system("mv Correlation.txt Correlation.ps Correlation.pdf Correlation.eps plots")
  call system("mv Correlation.gnu plots/Correlationgnus")

  endif !if plot at all
  end subroutine KL_Correlation








  subroutine KL_collect( nummatSegs,matLength,matType,j,&
                         numEigs,sigave,lamc,totLength,&
                         time,ntime )
  use genRealzvars, only: sig, lam, s, numRealz
  use KLvars,       only: gam, alpha, Ak, Eig, xi
  integer :: matType(:),j,numEigs,ntime,nummatSegs
  real(8) :: matLength(:),sigave,time(:),tt1,tt2
  real(8) :: lamc,totLength(2)

  integer :: i,k,curEig,doloop
  real(8) :: xitermtot,xl,xr,sigma,xiterm

  call cpu_time(tt1)

  do curEig=1,numEigs
    xl=0
    xiterm=0
    xitermtot=0

    do i=2,nummatSegs+1
      xl=matLength(i-1)   !set xl and xr for calculations
      xr=matLength(i)
      sigma=sig(matType(i-1))                    !set sig to the correct sig

      xiterm= (sigma-sigave)*&                 !actual calculation
              (lamc*sin(alpha(curEig)*xr)-cos(alpha(curEig)*xr)/alpha(curEig) &
              -lamc*sin(alpha(curEig)*xl)+cos(alpha(curEig)*xl)/alpha(curEig))

      xitermtot = xitermtot + xiterm
    enddo

    xi(j,curEig) = (Ak(curEig)/sqrt(Eig(curEig)))*xitermtot   !find resulting xi
  enddo


  do i=2,nummatSegs !collect total length data for Actual Co calculations
    totLength(matType(i-1))=totLength(matType(i-1))+matLength(i)-matLength(i-1)
  enddo

  call cpu_time(tt2)
  time(5) = time(5) + (tt2-tt1)

  end subroutine KL_collect









  subroutine KL_Cochart( numEigs,numSlice,P,sigave,lamc,&
                         avePath,totLength,pltCo,&
                         CoExp )
  !This subroutine calculates the ratio of the calculated variace (Co) using a chosen
  !number of eigenmodes to the total variance, which is equivalent to using all
  !eigenmodes.  The ratio will thus always be less than 1.  The close to 1 the ratio
  !is, the more efficient that approximation is.  This calculation is made using the
  !Co value given by input parameters (P(1) & P(2)), and by the generated actual
  !probabilities.  Since the eigenvalue contains the variance, this baseline Co
  !actually divides itself out, so that the efficiency by either method is the same.
  !This subroutine calculates both, then prints those that are chosen in the input.
  use genRealzvars, only: sig, s, numRealz
  use KLvars,       only: gam, alpha, Ak, Eig, pltCowhich, pltConumof
  integer :: numEigs,numSlice
  real(8) :: P(2),sigave,lamc,avePath(2),totLength(2),CoExp
  character(7) :: pltCo(4)

  integer :: curCS,curEig,twice,check
  real(8) :: slicesize,cumCo,sliceval(numSlice),CoEff(numEigs,numSlice)
  real(8) :: Co,CoAct,CoPerDiff,totLPer(2),tottotLength
  real(8) :: x,phi
  real(8),allocatable :: Coplotarray(:,:)

  allocate(Coplotarray(numSlice,pltConumof+1))
  Coplotarray = 0

  444 format("totLength(1):",f15.2,"    totLength(2):",f15.2,"    tottotLength:",f15.2)
  445 format("expLength(1):",f15.2,"    expLength(2):",f15.2,"    exptotLength:",f15.2)
  446 format("       CoAct:",f15.5)
  447 format("       CoExp:",f15.5,"       CoPerdiff:",f15.2," %")


  totLPer(1)=totLength(1)/(totLength(1)+totLength(2))
  totLPer(2)=totLength(2)/(totLength(1)+totLength(2))
  tottotLength=totLength(1)+totLength(2)
  write(*,*)
  write(*,444) totLength(1),totLength(2),tottotLength
  write(*,445) P(1)*s*numRealz,P(2)*s*numRealz,s*numRealz


  slicesize=s/numSlice     !prepare the places to slice and measure
  sliceval(1)=slicesize/2
  do curCS=2,numSlice
    sliceval(curCS)=sliceval(curCS-1)+slicesize
  enddo

  !CoExp     = P(1)       *P(2)       *(sig(1)-sig(2))**2 !already calced
  CoAct     = totLPer(1) *totLPer(2) *(sig(1)-sig(2))**2
  CoPerDiff = abs( ( CoExp - CoAct ) / CoExp ) * 100
  write(*,446) CoAct
  write(*,447) CoExp,CoPerDiff


  do twice=1,2

    if( twice==1 ) then
      Co=CoExp
    else
      Co=CoAct
      do curEig=1,numEigs
        Eig(curEig)=2*Co*lamc/(gam(curEig)**2+1)
      enddo
    endif

    do curCS=1,numSlice      !calculate Co(ours)/Co; for slices
      cumCo=0
      x=sliceval(curCS)
      do curEig=1,numEigs
        phi=Eigfunc(Ak(curEig),alpha(curEig),lamc,x)
        cumCo=cumCo+  Eig(curEig)*phi**2

        CoEff(curEig,curCS)=cumCo/Co
      enddo
    enddo


    440 format("#    slicenum    sliceval ")
    441 format("     CoEff",i5)
    442 format("     ",i5,"     ",f10.6)
    443 format("     ",f10.6)

    open(unit=14, file="CoEff.txt")
    write(14,440,advance="no")
    do curEig=1,numEigs
      write(14,441,advance="no") curEig
    enddo
    do curCS=1,numSlice       !print CoEff results to file
      write(14,*)
      write(14,442,advance="no") curCS,sliceval(curCS)
      do curEig=1,numEigs
        write(14,443,advance="no") CoEff(curEig,curCS)
      enddo
    enddo

    if( twice==1 ) then
      call system("mv CoEff.txt texts/CoEffExp.txt")
    else
      call system("mv CoEff.txt texts/CoEffAct.txt")
    endif

    do check=1,pltConumof  !load values to plot according to selected options
      if( pltCowhich(2,check)==twice .AND. pltCo(1) .NE. 'noplot' ) then
        curEig = pltCowhich(1,check)
        do curCS=1,numSlice
          Coplotarray(curCS,1) = sliceval(curCS)
          Coplotarray(curCS,check+1) = CoEff(curEig,curCS)
        enddo
      endif
    enddo

  enddo !end twice, loop for expected and actual Co values

  if(pltCo(1) .NE. 'noplot') then
    call generic_plotter( numSlice,pltConumof,Coplotarray,pltCo ) !plot

    call system("mv genericplot.txt plots/CoEffplot.txt")
    call system("mv genericplot.ps  plots/CoEffplot.ps")
    call system("mv genericplot.pdf plots/CoEffplot.pdf")
    call system("rm plots/CoEffplot.eps")
    call system("mv genericplot.eps plots/CoEffplot.eps")
  endif


  do curEig=1,numEigs  !return 'Eig' to original values
    Eig(curEig)=2*CoExp*lamc/(gam(curEig)**2+1)
  enddo

  end subroutine KL_Cochart











  subroutine KL_eval( binSmallBound,binLargeBound,&
                      numEigs,&
                      pltxiBins,pltxiBinsgauss,binSize,&
                      mostinBin )
  !This subroutine puts xi values in bins.  It plots for those chosen in the 
  !input file and also makes a .txt file containing PDFs of xi values for 
  !each Eigenvalue calculated.
  use genRealzvars, only: numRealz
  use KLvars, only: xi, binPDF, binBounds, pltxiBinswhich, pltxiBinsnumof, binNumof
  integer :: numEigs,mostinBin
  real(8) :: binSmallBound,binLargeBound
  real(8),allocatable :: binper(:,:)
  character(7) :: pltxiBins(4),pltxiBinsgauss

  integer :: k,j,i,binCounts(binNumof,numEigs+1),curEig,tnumRealz
  real(8) :: smallestxi,largestxi,xiOneD(numRealz)
  integer :: binCountsplotarray(binNumof,10),binCountsOneD(binNumof)
  real(8) :: probsum,meanxi,varxi,binSize
  real(8) :: pi=3.14159265358979d0
  real(8),allocatable :: binPDFplotarray(:,:)

  allocate(binPDF(binNumof,numEigs+1))
  allocate(binper(binNumof,numEigs+1))
  allocate(binPDFplotarray(binNumof,pltxiBinsnumof+1))
  allocate(binBounds(binNumof+1))

  binCounts = 0
  binBounds = 0
  mostinBin = 0

  !print xi to screen
  !do j=j,numRealz
  !  do curEig=1,numEigs
  !    480 format("xi( ):  ",f18.5)
  !    write(*,480) xi(j,curEig)
  !  enddo
  !enddo
  
  smallestxi=minval(xi)  !get largest and smallest xi values
  largestxi=maxval(xi)

  if( binSmallBound==0 ) then        !set bounds for bins if needed
    binSmallBound=floor(smallestxi)
  endif
  if( binLargeBound==0 ) then
    binLargeBound=ceiling(largestxi)
  endif
  482 format("  smallestxi:",f24.5,"        largestxi:",f24.5)
  write(*,*)
  write(*,*)
  write(*,482) smallestxi,largestxi
  write(*,*)
  binSize = ( binLargeBound - binSmallBound ) / ( binNumof - 1 )
  binCountsplotarray = 0
  binPDFplotarray = 0

  !This section prints (and plots) for Eigs chosen to plot
  if( pltxiBins(1) .NE. 'noplot' ) then  !plot using generic plotter
    do i=1,pltxiBinsnumof  !up to 4 options, collect data to plot
      curEig    = pltxiBinswhich(1,i)
      tnumRealz = pltxiBinswhich(2,i)

      do j=1,numRealz   !prepare xi and binCounts for binner
        xiOneD(j) = xi(j,curEig)
      enddo
      binCountsOneD = 0
      call store_in_bins( binSmallBound,binLargeBound,binNumof,&
                          binCountsOneD,binBounds,xiOneD,tnumRealz )
      do k=1,binNumof    !collect binCounts from binner
        binCountsplotarray(k,i+1) = binCountsOneD(k)
      enddo

      probsum = 0  !transition data to PDF form to plot
      do k=1,binNumof
        binPDFplotarray(k,1)  =binBounds(k)+binSize/2
        binPDFplotarray(k,i+1)=real(binCountsplotarray(k,i+1),8)/&
                               real(tnumRealz,8)/binSize
        if( pltxiBinsgauss .EQ. 'gauss' .AND. i==pltxiBinsnumof ) then
          binPDFplotarray(k,i+1)=exp(-binPDFplotarray(k,1)**2/2)/(sqrt(2*pi))
        endif
        probsum=probsum+binPDFplotarray(k,i+1)*binSize
      enddo
      481 format("  curEig:",i3,"  tnumRealz:",i11,"   probsum:",f6.3)
      write(*,481) curEig,pltxiBinswhich(2,i),probsum
    enddo

    call generic_plotter( binNumof,pltxiBinsnumof,binPDFplotarray,pltxiBins )

  !rename "generic" plotting files
    call system("mv genericplot.txt plots/xiBinsplot.txt")
    call system("mv genericplot.ps  plots/xiBinsplot.ps")
    call system("mv genericplot.pdf plots/xiBinsplot.pdf")
    call system("rm plots/xiBinsplot.eps")
    call system("mv genericplot.eps plots/xiBinsplot.eps")
  endif


  !This section prints for all Eigs, not only those to plot
  write(*,*)
  do curEig=1,numEigs
    do j=1,numRealz   !prepare xi and binCounts for binner
      xiOneD(j) = xi(j,curEig)
    enddo
    binCountsOneD = 0
    call store_in_bins( binSmallBound,binLargeBound,binNumof,&
                        binCountsOneD,binBounds,xiOneD,numRealz )
    do k=1,binNumof   !collect binCounts from binner, and track most
      binCounts(k,curEig+1) = binCountsOneD(k)
      if( binCountsOneD(k) > mostinBin) mostinBin = binCountsOneD(k)
    enddo

    probsum = 0
    meanxi  = 0
    varxi   = 0
    do k=1,binNumof
      binPDF(k,1) = binBounds(k)+binSize/2
      binPDF(k,curEig+1) = real(binCounts(k,curEig+1),8)/real(numRealz,8)/binSize
      binper(k,1) = binPDF(k,1)
      binper(k,curEig+1) = real(binCounts(k,curEig+1),8)/real(numRealz,8)*100.0d0
      probsum = probsum + binPDF(k,curEig+1) * binSize
      meanxi  = meanxi  + binPDF(k,curEig+1) * binSize * binPDF(k,1)
      varxi   = varxi   + binPDF(k,cureig+1) * binSize * binPDF(k,1) ** 2
    enddo
    483 format("  curEig:",i3,"   numRealz:",i11,"   probsum:",f6.3,"   meanxi:",f9.4,"   varxi:",f9.4)
    write(*,483) curEig,numRealz,probsum,meanxi,varxi
  enddo
  write(*,*)

  open(unit=1,file="xiBinsAll.txt")  !Print data to file
  open(unit=2,file="xiBinsper.txt")  !Print percentages to file
  do i=1,binNumof
    484 format(f12.6,"     ")
    485 format(f11.5," %    ")
    write(1,484,advance="no") binPDF(i,1)
    write(2,484,advance="no") binper(i,1)
    do j=2,numEigs+1
      write(1,484,advance="no") binPDF(i,j)
      write(2,485,advance="no") binper(i,j)
    enddo
    write(1,*)
    write(2,*)
  enddo
  close(unit=1)
  close(unit=2)
  call system("mv xiBinsAll.txt plots")
  call system("mv xiBinsper.txt plots")

  end subroutine KL_eval








  subroutine KL_Noise( numEigs,&
                       binSmallBound,binLargeBound,binSize,mostinBin,&
                       time,ntime )
  !This subroutine identifies the two largest peaks of a xi distribution,
  !then uses the variable "xi" to print which realization/eigenmode combos
  !yielded values outside of these two peaks... which were the "noise" under
  !investigation.  Data printed is Eigs for each are out, and how far from
  !the nearest "largest peak", neg for below bottom, pos for above top.
  use genRealzvars, only: numRealz
  use Klvars, only: xi, binPDF, binBounds, binNumof
  integer :: numEigs,mostinBin,ntime
  real(8) :: binSmallBound,binLargeBound
  real(8) :: binSize,time(:),tt1,tt2

  integer :: j,i,curEig,indofbin,w
  integer :: whererealz(numEigs,binNumof,mostinBin+1)
  real(8) :: peaks(numEigs,2,2) !eigenmode;peak 1 or 2;height of peak,loc of peak
  character(3) :: search

  call cpu_time(tt1)

  !This part tracks which realization was put in which bin, then prints
  !them all for each eigenmode to a text file
  !print *," mostinBin:",mostinBin
  whererealz = 0

  do curEig=1,numEigs !create array "whererealz", tracking each realization/eigenmode
    do j=1,numRealz
      do i=1,binNumof

        if( xi(j,curEig) > binBounds(i) .AND. xi(j,curEig) < binBounds(i+1) ) then
          !print *,"it's a fit!",xi(j,curEig),"between",binBounds(i),"and",binBounds(i+1)
          indofbin =   0
          search   = 'yes'
          do while ( search == 'yes' )
            indofbin = indofbin + 1
            if( whererealz(curEig,i,indofbin) == 0 ) then
              whererealz(curEig,i,indofbin) = j
              search = 'no'
            endif
          enddo
        endif

      enddo
    enddo
  !print *,
  enddo

  open(unit=14,file="whererealz.txt")  !Print data to file

  write(14,*) "begin printing results"  !print array "whererealz"
  do curEig=1,numEigs
    450 format("Eig ",i6,"  starting here")
    write(14,*)
    write(14,450) curEig
    do i=1,binNumof
      439 format("curEig;",i3,"  bin center;",f9.5,"  bin #;",i8,": ")
      write(14,439,advance='no') curEig,binPDF(i,1),i

      indofbin =   0
      search   = 'yes'
      do while ( search == 'yes' )
        indofbin = indofbin + 1
        if( whererealz(curEig,i,indofbin) /= 0 ) then
          430 format(i2,",")
          431 format(i3,",")
          432 format(i4,",")
          433 format(i5,",")
          434 format(i6,",")
          435 format(i7,",")
          436 format(i8,",")
          437 format(i9,",")
          438 format(i14,",")
          w = whererealz(curEig,i,indofbin)
          if(w<10                         ) write(14,430,advance='no') w
          if(w<100       .AND. w>=10      ) write(14,431,advance='no') w
          if(w<1000      .AND. w>=100     ) write(14,432,advance='no') w
          if(w<10000     .AND. w>=1000    ) write(14,433,advance='no') w
          if(w<100000    .AND. w>=10000   ) write(14,434,advance='no') w
          if(w<1000000   .AND. w>=100000  ) write(14,435,advance='no') w
          if(w<10000000  .AND. w>=1000000 ) write(14,436,advance='no') w
          if(w<100000000 .AND. w>=10000000) write(14,437,advance='no') w
          if(w>100000000                 ) write(14,438,advance='no') w
        else
          search = 'no'
        endif
      enddo

      write(14,*)
    enddo
    write(14,*)
  enddo

  close(unit=14)
  call system("mv whererealz.txt plots")

  call cpu_time(tt2)
  time(4) = time(4) + (tt2-tt1)

!  allocate(xi(numRealz,numEigs))














!  peaks = 0

  !Identify two largest peaks for each eigenmode
!  do curEig=1,numEigs
!    do i=1,binNumof
!      if( binPDF(i,curEig+1) > peaks(curEig,1,1) ) then
!        peaks(curEig,2,1) = peaks(curEig,1,1) !top to second
!        peaks(curEig,2,2) = peaks(curEig,1,2)
!        peaks(curEig,1,1) = binPDF(i,curEig+1) !new top
!        peaks(curEig,1,2) = binPDF(i,1)
!      elseif( binPDF(i,curEig+1) > peaks(curEig,2,1) ) then
!        peaks(curEig,2,1) = binPDF(i,curEig+1) !new second
!        peaks(curEig,2,2) = binPDF(i,1)
!      endif
!    enddo
!  enddo

!print *,"height",peaks(1,1,1),"loc",peaks(1,1,2),"height",peaks(1,2,1),"loc",peaks(1,2,2)
!print *,"height",peaks(2,1,1),"loc",peaks(2,1,2),"height",peaks(2,2,1),"loc",peaks(2,2,2)
!print *,"height",peaks(3,1,1),"loc",peaks(3,1,2),"height",peaks(3,2,1),"loc",peaks(3,2,2)
!print *,"height",peaks(4,1,1),"loc",peaks(4,1,2),"height",peaks(4,2,1),"loc",peaks(4,2,2)
!print *,"height",peaks(5,1,1),"loc",peaks(5,1,2),"height",peaks(5,2,1),"loc",peaks(5,2,2)

!  allocate(binPDF(binNumof,numEigs+1))

  end subroutine KL_Noise





end module KLresearch

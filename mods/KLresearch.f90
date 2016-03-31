module KLresearch
  use utilities
  implicit none

CONTAINS
  ! print statements in this module use # 400-499


  subroutine KL_eigenvaluemain()
  !This subroutine is the interface for 'KL_eigenvalue' so that values can be passed to
  !it by reference.
  use KLvars, only: alphas1,Aks1,Eigs1,numEigss1, alphaa1,Aka1,Eiga1,numEigsa1,&
                    alphas2,Aks2,Eigs2,numEigss2, alphaa2,Aka2,Eiga2,numEigsa2,&
                    fls1, fla1, fls2, fla2
  use genRealzvars, only: lamcs1,lamca1,lamcs2,lamca2

  if(fls1) call KL_eigenvalue(alphas1,Aks1,Eigs1,numEigss1,lamcs1)
  if(fla1) call KL_eigenvalue(alphaa1,Aka1,Eiga1,numEigsa1,lamca1)
  if(fls2) call KL_eigenvalue(alphas2,Aks2,Eigs2,numEigss2,lamcs2)
  if(fla2) call KL_eigenvalue(alphaa2,Aka2,Eiga2,numEigsa2,lamca2)

  end subroutine KL_eigenvaluemain


  function Eigenvalue( gam,lamc )
  !Solves value of eigenvalues based on transcendental solution gamma.
  real(8) :: gam, lamc, Eigenvalue
  Eigenvalue = 2d0 * lamc / (gam**2 + 1d0)
  end function Eigenvalue



  subroutine KL_eigenvalue(alpha,Ak,Eig,numEigs,lamc)
  !This subroutine: 1) calculates some initial values used here
  !2) Solves the transcendental equation which yields gamma
  !3) From gamma solves: alpha, lambda (Eigenvalue), & the normalization const A_k
  !4) Prints and plots Eigenfunctions if input specifies
  !5) Calculates the percent of mean standard error maintained
  use genRealzvars, only: s, aves1, avea1, GBaves1, GBavea1, chgeomtype, GBs, GBvars1, GBvara1
  use KLvars,       only: levsrefEig, pltEigf, pltEigfwhich, pltEigfnumof, numSlice, &
                          lamctypes1,cheftypes1,numNystroms1,eigvecss1,Eigs1, &
                          lamctypea1,cheftypea1,numNystroma1,eigvecsa1,Eiga1
  use KLconstruct, only: Eigfunc

  real(8) :: alpha(:),Ak(:),Eig(:),lamc
  integer :: numEigs

  real(8) :: stepGam=0 !if 0 code chooses
  integer :: l,level,curEig,i,j
  real(8) :: refstepGam,TT,curGam,tc
  real(8) :: absdiff,absdiff_1=0,absdiff_2=0,testval(numEigs),sliceSize
  real(8),allocatable :: Eigfplotarray(:,:),gam(:)

  TT = s/lamc
  if( stepGam==0 ) stepGam  =1/TT/50
  allocate(gam(numEigs))
  gam = 0.0d0
  allocate(Eigfplotarray(numSlice,numSlice+1))

  !Initial guesses for Gam
  420 format(i7,"     ",f15.9," ",f19.14)

  write(*,*) " Calculating Eigenvalues for computation and varaince "
  write(*,*) "          Eigindx    Gam vals       Eig vals         tol check         tol "

  !solve gamma approximate values
  curGam = 0d0
  curEig = 1
  open(unit=11, file="absdiffGam.txt")
  do
    absdiff=abs( tan(curGam*TT)-(2d0*curGam)/(curGam**2-1d0) )
    write(11,420) curEig,curGam,absdiff
    if( absdiff_2>absdiff_1 .and. absdiff>absdiff_1 ) then
      !write(*,420) curEig,curGam,absdiff_1
      gam(curEig)=curGam
      if(curEig==numEigs) exit
      curEig=curEig+1
    endif
    absdiff_2=absdiff_1
    absdiff_1=absdiff
    curGam=curGam+stepGam
  enddo
  close(unit=11)
  call system("mv absdiffGam.txt texts/")

  !refine gamma values and solve eigenvalues
  do curEig=1,numEigs                  !loop each Eigenvalue
    refstepGam=stepGam

    do level=1,levsrefEig                 !loop levels of refinement
      curGam=gam(curEig)-11d0*refstepGam
      refstepGam=refstepGam/10
      absdiff_1=0
      absdiff_2=0

      do l=1,220                            !loop through range & test
        absdiff=abs( tan(curGam*TT)-(2d0*curGam)/(curGam**2-1d0) )
        if( absdiff_2>absdiff_1 .AND. absdiff>absdiff_1 ) then
          !write(*,421) curEig,curGam,absdiff_1,level
          !if(level==levsrefEig) write(*,421) curEig,curGam,absdiff_1,level
          gam(curEig)=curGam
        endif
        absdiff_2=absdiff_1
        absdiff_1=absdiff
        curGam=curGam+refstepGam
      enddo
    enddo
    Eig(curEig) = Eigenvalue( gam(curEig),lamc )
  enddo


  write(*,*) "     Int test"
  do curEig=1,numEigs
    !Calc other values like alpha, norm const (Ak), eigenvalue, etc.
    alpha(curEig)   =gam(curEig)/lamc
    Ak(curEig)      =sqrt(1d0/(  s/2d0*(gam(curEig)**2+1d0)+lamc  ))
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


  write(*,*)
  write(*,*) "    Eigenvalues, their contributions, and KL maintained variance"
  write(*,*) "  Eigindx       Eig vals        sqrt(Eig)     [ differential / cumulative ] % mean sqr err maint"
  426 format(i7,"     ",f15.9,"  ",f13.9,"     ",f13.9,"  ",f13.9)
!  if(chgeomtype=='binary') then
!    tc = aves1/(aves1+avea1)
!  elseif(chgeomtype=='contin') then
!    tc = GBaves1/(GBaves1+GBavea1)
!  endif
!  do curEig=1,numEigs
!    write(*,426) curEig,Eig(curEig),sqrt(Eig(curEig)),&
!                (   tc   *merge(Eig(curEig),0d0,curEig<=numEigss1) + &
!                 (1d0-tc)*merge(Eig(curEig),0d0,curEig<=numEigsa1)) / s,&
!                (   tc   *sum( Eig(:min(curEig,numEigss1)) ) + &
!                 (1d0-tc)*sum( Eig(:min(curEig,numEigsa1)) )) /s
!  enddo

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
        Eigfplotarray(i,j+1) = Eigfunc( Ak(curEig),alpha(curEig),lamc,Eigfplotarray(i,1),0,        &
                                        lamctypes1,cheftypes1,numNystroms1,GBs,eigvecss1(curEig,:),&
                                        GBaves1,GBvars1,Eigs1(curEig)  )

      enddo
    enddo

    call generic_plotter( numSlice,pltEigfnumof,Eigfplotarray,pltEigf )

    call system("mv genericplot.txt plots/Eigfplot/Eigfplot.txt")
    call system("mv genericplot.ps  plots/Eigfplot/Eigfplot.ps")
    call system("mv genericplot.pdf plots/Eigfplot/Eigfplot.pdf")
  endif


  end subroutine KL_eigenvalue




  subroutine KL_collect
  !This subroutine calculates xi values from binary Markov realizations.  It does 
  !this by integrating over a variable related only to geometry-based material
  !properties so that it can later be used for any cross section material property.
  use rngvars, only: rngappnum, rngstride
  use genRealzvars, only: sig, numRealz, nummatSegs, lamcs1, matType, matLength, P
  use KLvars, only: alphas1, Aks1, Eigs1, xi, numEigss1
  use MCvars, only: trannprt
  use genRealz, only: genbinaryReal
  use timeman, only: initialize_t1, timeupdate
  use mcnp_random, only: RN_init_particle
  integer :: i,j,curEig,lsig,ssig
  real(8) :: xitermtot,xl,xr,xiterm,  hilowterm,aveterm
  character(5) :: flKLtype = 'KLcol'

  call initialize_t1

  write(*,*) "Starting method: ",flKLtype  
  lsig = merge(1,2,sig(1)>sig(2))
  ssig = merge(1,2,sig(1)<sig(2))
  !advance rng
  call RN_init_particle( int(rngappnum*rngstride,8) )
  rngappnum = rngappnum + 1

  do j=1,numRealz
    call genbinaryReal( j )
    do curEig=1,numEigss1
      xl       =0d0
      xiterm   =0d0
      xitermtot=0d0

      do i=2,nummatSegs+1
        xl=matLength(i-1)                                         !set xl and xr for calculations
        xr=matLength(i)
        hilowterm = merge(sqrt(P(ssig)/P(lsig)),-sqrt(P(lsig)/P(ssig)),matType(i-1)==lsig)
        aveterm   = 0d0

        xiterm= (hilowterm-aveterm)*&                                  !actual calculation
                (lamcs1*sin(alphas1(curEig)*xr)-cos(alphas1(curEig)*xr)/alphas1(curEig) &
                -lamcs1*sin(alphas1(curEig)*xl)+cos(alphas1(curEig)*xl)/alphas1(curEig))

        xitermtot = xitermtot + xiterm
      enddo
      xi(j,curEig) = (Aks1(curEig)/sqrt(Eigs1(curEig)))*xitermtot     !find resulting xi
    enddo

    if(mod(j,trannprt)==0) call timeupdate( 'KL_collect',j,numRealz )

  enddo

  end subroutine KL_collect





  subroutine KL_binrandvarvals
  !This subroutine puts xi values in bins.  It plots for those chosen in the 
  !input file and also makes a .txt file containing PDFs of xi values for 
  !each Eigenvalue calculated.
  use genRealzvars, only: numRealz
  use KLvars, only: xi, binPDF, binBounds, pltxiBinswhich, pltxiBinsnumof, binNumof, &
                    numEigss1, mostinBin, binSmallBound, binLargeBound, binSize, pltxiBins, &
                    pltxiBinsgauss
  real(8),allocatable :: binper(:,:)

  integer :: k,j,i,binCounts(binNumof,numEigss1+1),curEig,tnumRealz
  real(8) :: smallestxi,largestxi,xiOneD(numRealz)
  integer :: binCountsplotarray(binNumof,10),binCountsOneD(binNumof)
  real(8) :: probsum,meanxi,varxi
  real(8) :: pi=3.14159265358979d0
  real(8),allocatable :: binPDFplotarray(:,:)

  allocate(binPDF(binNumof,numEigss1+1))
  allocate(binper(binNumof,numEigss1+1))
  allocate(binPDFplotarray(binNumof,pltxiBinsnumof+1))
  allocate(binBounds(binNumof+1))

  binCounts = 0
  binBounds = 0
  mostinBin = 0

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

    call system("mv genericplot.txt plots/xiBinsplot/xiBinsplot.txt")
    call system("mv genericplot.ps  plots/xiBinsplot/xiBinsplot.ps")
    call system("mv genericplot.pdf plots/xiBinsplot/xiBinsplot.pdf")
  endif


  !This section prints for all Eigs, not only those to plot
  write(*,*)
  do curEig=1,numEigss1
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
    do j=2,numEigss1+1
      write(1,484,advance="no") binPDF(i,j)
      write(2,485,advance="no") binper(i,j)
    enddo
    write(1,*)
    write(2,*)
  enddo
  close(unit=1)
  close(unit=2)
  call system("mv xiBinsAll.txt plots/xiBinsplot")
  call system("mv xiBinsper.txt plots/xiBinsplot")

  end subroutine KL_binrandvarvals




  subroutine KL_Correlation
  !This subroutine calculates Correlation between two points in a
  !realization based upon the expected value, and the observed 
  !value (function of Eigenfunctions and values).
  !It then plots in 3D if user has specified.
  use genRealzvars, only: s, lamcs1, aves1, avea1, GBaves1, GBavea1, chgeomtype, GBs, GBvara1
  use KLvars, only: alphas1, Aks1, Eigs1, Corrnumpoints, Corropts, numEigss1, numEigsa1, & 
                    lamctypea1,cheftypea1,numNystroma1,eigvecsa1,Eiga1
  use KLconstruct, only: Eigfunc

  integer :: x,y,curEig
  real(8) :: stepsize,curx,cury,Eigfx,Eigfy,tc
  real(8) :: Correxpect(Corrnumpoints,Corrnumpoints)
  real(8) :: Corryield(Corrnumpoints,Corrnumpoints)
  character(17) :: Corrchar = '"Correlation.txt"' 
  character(22) :: Corre    = '"Correlation Expected"'
  character(22) :: Corry    = '"Correlation Yielded "'

  Correxpect = 0  !initialize
  Corryield  = 0

  stepsize = s/(Corrnumpoints)  !set up stepsize

  if(chgeomtype=='binary') then
    tc = aves1/(aves1+avea1)
  elseif(chgeomtype=='contin') then
    tc = GBaves1/(GBaves1+GBavea1)
  endif

  do x=1,Corrnumpoints  !cycle through x and y
    if(x==1) curx=stepsize/2 
    do y=1,Corrnumpoints
      if(y==1) cury=stepsize/2

      !calc Correxpect
      Correxpect(x,y) = exp( - abs(curx-cury)/lamcs1 )

      !sum Eigs to calc Corryield
      do curEig=1,numEigss1
        Eigfx = Eigfunc( Aks1(curEig),alphas1(curEig),lamcs1,curx,0,                &
                         lamctypea1,cheftypea1,numNystroma1,GBs,eigvecsa1(curEig,:),&
                         GBavea1,GBvara1,Eiga1(curEig)  )
        Eigfy = Eigfunc( Aks1(curEig),alphas1(curEig),lamcs1,cury,0,                &
                         lamctypea1,cheftypea1,numNystroma1,GBs,eigvecsa1(curEig,:),&
                         GBavea1,GBvara1,Eiga1(curEig)  )
        Corryield(x,y) = Corryield(x,y) + &
                         merge(    tc   *(Eigs1(curEig) * Eigfx * Eigfy),0d0,curEig<=numEigss1 ) + &
                         merge( (1d0-tc)*(Eigs1(curEig) * Eigfx * Eigfy),0d0,curEig<=numEigsa1 )
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
    call system("cat Correlationfirst.gnu plots/Correlation/Corrps_pause.txt > Correlation.gnu")
  else
    call system("cat Correlationfirst.gnu plots/Correlation/Corrps.txt > Correlation.gnu")
  endif

  call system("gnuplot Correlation.gnu") !plot using gnuplot

  call system("ps2pdf Correlation.ps")         !convert to pdf

  call system("rm Correlationfirst.gnu")                       !clean up space
  call system("mv Correlation.txt Correlation.ps Correlation.pdf plots/Correlation")
  call system("mv Correlation.gnu plots/Correlation")
  end subroutine KL_Correlation




  subroutine KL_Cochart
  !This subroutine calculates the variance normalized to 1 at each point in the domain.
  !The closer to 1 the ratio is, the more efficient that approximation is.
  use genRealzvars, only: s, numRealz, P, lamcs1, totLength, chgeomtype, aves1, avea1, &
                          GBaves1, GBavea1, GBs, GBvara1
  use KLvars,       only: alphas1, Aks1, Eigs1, pltCowhich, pltConumof, numSlice, &
                          pltCo, numEigss1, numEigsa1, lamctypea1,cheftypea1,numNystroma1,eigvecsa1,Eiga1
  use KLconstruct, only: Eigfunc

  integer :: curCS,curEig,check
  real(8) :: slicesize,cumCo,sliceval(numSlice),CoEff(numEigss1,numSlice)
  real(8) :: totLPer(2),tottotLength
  real(8) :: x,phi,tc
  real(8),allocatable :: Coplotarray(:,:)

  allocate(Coplotarray(numSlice,pltConumof+1))
  Coplotarray = 0

  !This is a statistical test of the reconstructions
  444 format("totLength(1):",f15.2,"    totLength(2):",f15.2,"    tottotLength:",f15.2)
  445 format("expLength(1):",f15.2,"    expLength(2):",f15.2,"    exptotLength:",f15.2)

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

  if(chgeomtype=='binary') then
    tc = aves1/(aves1+avea1)
  elseif(chgeomtype=='contin') then
    tc = GBaves1/(GBaves1+GBavea1)
  endif

  do curCS=1,numSlice      !calculate Co(ours)/Co; for slices
    cumCo=0
    x=sliceval(curCS)
    do curEig=1,numEigss1
      phi=Eigfunc( Aks1(curEig),alphas1(curEig),lamcs1,x,0,                   &
                   lamctypea1,cheftypea1,numNystroma1,GBs,eigvecsa1(curEig,:),&
                   GBavea1,GBvara1,Eiga1(curEig)  )
      cumCo=cumCo+ merge(   tc   *Eigs1(curEig)*phi**2,0d0,curEig<=numEigsa1) &
                 + merge((1d0-tc)*Eigs1(curEig)*phi**2,0d0,curEig<=numEigss1)
      CoEff(curEig,curCS)=cumCo
    enddo
  enddo


  440 format("#    slicenum    sliceval ")
  441 format("     CoEff",i5)
  442 format("     ",i5,"     ",f10.6)
  443 format("     ",f10.6)

  open(unit=14, file="CoEff.txt")
  write(14,440,advance="no")
  do curEig=1,numEigss1
    write(14,441,advance="no") curEig
  enddo
  do curCS=1,numSlice       !print CoEff results to file
    write(14,*)
    write(14,442,advance="no") curCS,sliceval(curCS)
    do curEig=1,numEigss1
      write(14,443,advance="no") CoEff(curEig,curCS)
    enddo
  enddo

  call system("mv CoEff.txt texts/CoEffExp.txt")

  do check=1,pltConumof  !load values to plot according to selected options
    if( pltCo(1) .NE. 'noplot' ) then
      curEig = pltCowhich(check)
      do curCS=1,numSlice
        Coplotarray(curCS,1) = sliceval(curCS)
        Coplotarray(curCS,check+1) = CoEff(curEig,curCS)
      enddo
    endif
  enddo

  call generic_plotter( numSlice,pltConumof,Coplotarray,pltCo ) !plot

  call system("mv genericplot.txt plots/CoEffplot/CoEffplot.txt")
  call system("mv genericplot.ps  plots/CoEffplot/CoEffplot.ps")
  call system("mv genericplot.pdf plots/CoEffplot/CoEffplot.pdf")

  end subroutine KL_Cochart




end module KLresearch

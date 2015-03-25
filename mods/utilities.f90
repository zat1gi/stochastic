module utilities
  implicit none
CONTAINS
  ! print statements in this module use # 1000-1100

  subroutine calc_time_p( t1,t2,runtime )
  !Calculates runtime in minutes, does print
  real(8) :: t1,t2,runtime
  1000 format("  total runtime:",f7.2," min")

  call cpu_time(t2)
  runtime=(t2-t1)/60
  write(*,1000) runtime
  end subroutine calc_time_p
 



  subroutine calc_time( t1,t2,runtime )
  !Calculates runtime in minutes, does not print
  real(8) :: t1,t2,runtime

  call cpu_time(t2)
  runtime=(t2-t1)/60.0d0
  end subroutine calc_time




  subroutine store_in_bins( binSmallBound,binLargeBound,binNumof,&
             binCounts,binBounds,vararray,numVar )
  !Bins values between binSmallBound and binLargeBound with binNumof bins
  !vararray is what you are binning, with numVar as the length of vararray
  !binBounds is bounds with binNumof+1 size, counts in binCounts with binNumof as size
  integer, intent(in)  :: binNumof,numVar
  integer, intent(out) :: binCounts(binNumof)
  real(8), intent(in)  :: binSmallBound,binLargeBound,vararray(numVar)
  real(8), intent(out) :: binBounds(binNumof+1)

  integer :: binToStore,j,k
  real(8) :: binSize

  !set up bins
  binSize = ( binLargeBound - binSmallBound ) / ( binNumof - 1 )
  binBounds(1)=binSmallBound
  do k=2,binNumof+1
    binBounds(k) = binBounds(k-1) + binSize
  enddo

  !store in bins
  do j=1,numVar
    binToStore = ceiling(( vararray(j) - binSmallBound ) / binSize )
    binCounts(binToStore) = binCounts(binToStore) + 1
  enddo

  end subroutine store_in_bins





  subroutine generic_plotter( lengthCols,numCols,plotarray,plotoptions )
  !generic_plotter makes plots with the only nessesary imputs of the number
  !of columns, the length of columns, and an array with x-variables in the 
  !first column, and y-variables of up to 4 sets in the subsequent columns
  !The array must be (actualcolumn,columnnumber), and plotoptions must be
  !included.  See input file for clarification on plotoptions
  integer :: lengthCols,numCols
  real(8),allocatable :: plotarray(:,:)
  character(7) :: plotoptions(4)

  integer :: i,j

  open(unit=1,file="genericplot.txt")

  1012 format("              ",A7)
  write(1,1012,advance="no") plotoptions(3)
  do i=1,numCols
    1013 format("                 ",A7)
    write(1,1013,advance="no") plotoptions(4)
  enddo
  write(1,*)

  do i=1,lengthCols
    1014 format(f20.6,"     ")
    write(1,1014,advance="no") plotarray(i,1)
    do j=2,numCols+1
      write(1,1014,advance="no") plotarray(i,j)
    enddo
    write(1,*)
  enddo

  close(unit=1)


  if( plotoptions(1)=='plot' ) then  !choose the right generic plotter and plot
    if( plotoptions(2)=='lines' ) then
      if(numCols==1) call system("gnuplot plots/plottinggnus/generic.1.lines.gnu")
      if(numCols==2) call system("gnuplot plots/plottinggnus/generic.2.lines.gnu")
      if(numCols==3) call system("gnuplot plots/plottinggnus/generic.3.lines.gnu")
      if(numCols==4) call system("gnuplot plots/plottinggnus/generic.4.lines.gnu")
    elseif( plotoptions(2)=='points' ) then
      if(numCols==1) call system("gnuplot plots/plottinggnus/generic.1.points.gnu")
      if(numCols==2) call system("gnuplot plots/plottinggnus/generic.2.points.gnu")
      if(numCols==3) call system("gnuplot plots/plottinggnus/generic.3.points.gnu")
      if(numCols==4) call system("gnuplot plots/plottinggnus/generic.4.points.gnu")
    elseif( plotoptions(2)=='hist' ) then
      if(numCols==1) call system("gnuplot plots/plottinggnus/generic.1.hist.gnu")
      if(numCols==2) call system("gnuplot plots/plottinggnus/generic.2.hist.gnu")
      if(numCols==3) call system("gnuplot plots/plottinggnus/generic.3.hist.gnu")
      if(numCols==4) call system("gnuplot plots/plottinggnus/generic.4.hist.gnu")
    endif
  endif

  if( plotoptions(1)=='preview' ) then
    if( plotoptions(2)=='lines' ) then
      if(numCols==1) call system("gnuplot plots/plottinggnus/generic.1p.lines.gnu")
      if(numCols==2) call system("gnuplot plots/plottinggnus/generic.2p.lines.gnu")
      if(numCols==3) call system("gnuplot plots/plottinggnus/generic.3p.lines.gnu")
      if(numCols==4) call system("gnuplot plots/plottinggnus/generic.4p.lines.gnu")
    elseif( plotoptions(2)=='points' ) then
      if(numCols==1) call system("gnuplot plots/plottinggnus/generic.1p.points.gnu")
      if(numCols==2) call system("gnuplot plots/plottinggnus/generic.2p.points.gnu")
      if(numCols==3) call system("gnuplot plots/plottinggnus/generic.3p.points.gnu")
      if(numCols==4) call system("gnuplot plots/plottinggnus/generic.4p.points.gnu")
    elseif( plotoptions(2)=='hist' ) then
      if(numCols==1) call system("gnuplot plots/plottinggnus/generic.1p.hist.gnu")
      if(numCols==2) call system("gnuplot plots/plottinggnus/generic.2p.hist.gnu")
      if(numCols==3) call system("gnuplot plots/plottinggnus/generic.3p.hist.gnu")
      if(numCols==4) call system("gnuplot plots/plottinggnus/generic.4p.hist.gnu")
    endif
  endif

  call system("ps2pdf genericplot.ps genericplot.pdf")  !convert to pdf

  end subroutine generic_plotter







  subroutine select_from_PDF( binnedPDF,numbins,curcat,&
                              chosenvalue,rand )
  !This subroutine accepts a random number, then uses that number to 
  !sample from a PDF of bins.  It is currently set up where the first "category"
  !is the value represented by that bin, and each other category is a collection
  !of probabilities for choosing from the first, thus the size numcategories+1
  !numcategories+1 has been dynamically allocated, make sure it's correct!
  !This subroutine assumes that each bin is the same size
  integer :: numbins,curcat
  real(8) :: chosenvalue,rand
  real(8),allocatable :: binnedPDF(:,:)

  integer :: j
  real(8) :: binsize,newpercount,oldpercount

  binsize     = abs( binnedPDF(2,1) - binnedPDF(1,1) )
  newpercount = 0
  oldpercount = 0

  do j=1,numbins
    newpercount = oldpercount + binnedPDF(j,curcat+1) * binsize
    if( oldpercount<rand .AND. newpercount>=rand ) then
      chosenvalue = binnedPDF(j,1)
    endif
    oldpercount = newpercount
  enddo

  end subroutine select_from_PDF





  subroutine mean_and_var_p( values,numofvalues,mean,var )
  !This subroutine calculates mean and variance for a population, not a sample of
  !the population
  real(8) :: values(:),mean,var
  integer :: numofvalues

  integer :: i

  mean=0
  do i=1,numofvalues
    mean=mean+values(i)
  enddo
  mean=mean/numofvalues

  var=0
  do i=1,numofvalues
    var= var + (mean-values(i))**2
  enddo
  var=var/numofvalues

  end subroutine mean_and_var_p



  subroutine mean_and_var_s( values,numofvalues,mean,var )
  !This subroutine calculates mean and variance for a sample from a population,
  !not the whole population
  real(8) :: values(:),mean,var
  integer :: numofvalues

  integer :: i

  mean=0
  do i=1,numofvalues
    mean=mean+values(i)
  enddo
  mean=mean/numofvalues

  var=0
  do i=1,numofvalues
    var= var + (mean-values(i))**2
  enddo
  var=var/(numofvalues-1)

  end subroutine mean_and_var_s



  function mean(values,numofvalues)
  real(8) :: values(:),mean
  integer :: numofvalues,i

  mean=0
  do i=1,numofvalues
    mean=mean+values(i)
  enddo
  mean=mean/numofvalues

  end function mean



  function var_s(values,numofvalues)
  !This function calculates the variance for a sample from a population,
  !not the whole population
  real(8) :: values(:),mean,var_s
  integer :: numofvalues

  integer :: i

  mean=0
  do i=1,numofvalues
    mean=mean+values(i)
  enddo
  mean=mean/numofvalues

  var_s=0
  do i=1,numofvalues
    var_s= var_s + (mean-values(i))**2
  enddo
  var_s=var_s/numofvalues

  end function var_s




  function var_p(values,numofvalues)
  !This function calculates the variance for a whole population,
  !not a sample of the whole population
  real(8) :: values(:),mean,var_p
  integer :: numofvalues

  integer :: i

  mean=0
  do i=1,numofvalues
    mean=mean+values(i)
  enddo
  mean=mean/numofvalues

  var_p=0
  do i=1,numofvalues
    var_p= var_p + (mean-values(i))**2
  enddo
  var_p=var_p/(numofvalues-1)

  end function var_p



  function dev_s(values,numofvalues)
  !Standard deviation for a sample from a population
  real(8) :: values(:),var,dev_s
  integer :: numofvalues

  var = var_s(values,numofvalues)
  dev_s = sqrt(var)
  end function dev_s


  function dev_p(values,numofvalues)
  !Standard deviation for a whole population
  real(8) :: values(:),var,dev_p
  integer :: numofvalues

  var = var_p(values,numofvalues)
  dev_p = sqrt(var)
  end function dev_p




  function erfi(z)
  !Solves value for inverse error function to certain tolerance using
  !Maclaurin series.  Alerts user if not solved to tolerance due to 
  !numerical limitations.
  integer :: k,m
  real(8) :: z,erfi,erfinew
  real(8) :: tol = 0.00000000000000001d0
  real(8) :: pi  = 3.14159265358979323846d0
  real(8) :: largest = 1.6933d308
  integer :: most    = 2999
  logical :: fliterate
  real(8), dimension(:), allocatable :: c

  !screen valid function input
  if(z<=-1d0 .or. z>=1d0) stop "Input must be in interval (-1,1)"

  !initialize
  fliterate = .true.
  if(z==0d0) then
    erfi = 0d0
    fliterate = .false.
  endif
  allocate(c(-1:most))
  c     = 0d0
  c(-1) = 1d0
  erfi    = 0d0
  erfinew = 0d0

  k=0
  do while(fliterate)
    !solve needed coefficient
    c(k) = 0d0
    if(k/=0) then
      do m=0,k-1
        c(k) = c(k) + c(m)*c(k-1-m)/((real(m,8)+1d0)*(2d0*real(m,8)+1d0))
      enddo
    elseif(k==0) then
      c(k) = 1d0
    endif
    !use coefficient to solve next term
    erfinew = c(k)/(2d0*real(k,8)+1d0)*(sqrt(pi)/2d0*z)**(2d0*real(k,8)+1d0)
    !add term and/or stop loop
    if(erfinew>largest .or. erfinew<-largest .or. k>most) then
      fliterate = .false.
      print *,"erfi not solved below tolerance due to numerics"
    else
      if(erfinew/erfi<tol .and. k>0) fliterate = .false.
      erfi = erfi + erfinew
    endif
    k=k+1
  enddo
  deallocate(c)
  end function erfi



  subroutine Thomas_alg( n,a,b,c,d,x )
  !The thomas algorithm solves a tri-diagonal matrix.
  !This is the crux of a 1D diffusion model.
  integer               :: n
  real(8), dimension(n) :: a,b,c,d,x

  integer               :: i,j,checkalg

  checkalg=2
  if( checkalg==1 ) then  !a test for the Thomas algorithm
    a=(/   0,  -1,  -1,  -1/)
    b=(/2.04,2.04,2.04,2.04/)
    c=(/  -1,  -1,  -1,   0/)
    x=0
    d=(/40.8, 0.8, 0.8, 200.8/)
  endif !end of loading test values

  c(1)=    c(1)  /  b(1)
  d(1)=    d(1)  /  b(1)
  do i=2,n
    c(i)=  c(i)                    /  (  b(i)-c(i-1)*a(i)  )
    d(i)=  (  d(i)-d(i-1)*a(i)  )  /  (  b(i)-c(i-1)*a(i)  )
  enddo

  x(n)=d(n)
  do j=1,n-1
    i=n-j
    x(i)=  d(i)  -  c(i)*x(i+1)
  enddo

  if( checkalg==1 ) then !load and print correct values for input
    100 format("x(",i3,"):",f11.5)
    110 format("x(",i3,"):",f9.3)
    write(*,*)
    do i=1,n
      write(*,100) i,x(i)
    enddo
    write(*,*)
    x=(/65.970,93.778,124.538,159.48/)
    write(*,*) "Values (to 3 places) should be:"
    do i=1,n
      write(*,110) i,x(i)
    enddo
    write(*,*)
  endif !end of printing expected values for test

  end subroutine Thomas_alg




  subroutine gauss_leg_quad( numangs,nodes,wgts,leftb,rightb )
  !Solves nodes and weights with gauss-legendre quadrature, using executable.
  !Executable "gauss_leg_quad.exe", ensure in correct location.
  integer :: numangs
  real(8) :: leftb,rightb
  real(8),allocatable :: nodes(:),wgts(:)

  integer :: i

  !Allocate nodes and wgts
  allocate(nodes(numangs))
  allocate(wgts(numangs))
  nodes = 0.0
  wgts  = 0.0

  !Create driver script for executable
  open(unit=3,file="gl_run.sh")
  102 format("./auxiliary/quad/gauss_leg_quad.out ",i6," ",f12.6," ",f12.6," gauss_leg")
  write(3,102) numangs,leftb,rightb
  close(unit=3)

  !Use driver to run executable
  call system("chmod u+x gl_run.sh")
  call system("./gl_run.sh > gl_run.out")
  call system("rm gl_run.out")

  !Get data from output of executable
  open(unit=4,file="gauss_leg_x.txt")
  open(unit=5,file="gauss_leg_w.txt")
  do i=1,numangs
    read(4,*) nodes(i)
    read(5,*) wgts(i)
  enddo
  close(unit=4)
  close(unit=5)

  !Put all files into folder "quad"
  call system("mv gl_run.sh gauss_leg_r.txt gauss_leg_x.txt gauss_leg_w.txt auxiliary/quad")

  end subroutine gauss_leg_quad





end module utilities

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

  binCounts = 0
  binBounds = 0d0

  !set up bins
  binSize = ( binLargeBound - binSmallBound ) / ( binNumof )
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



  subroutine mean_and_var_wgt( wgts,values,mean,var )
  !This subroutine calculates the mean and variance for a data set with corresponding
  !weights.  The original purpose for this is use in stochastic collocation 
  !quantity of interest evaluation.
  real(8),intent(in) :: wgts(:),values(:)
  real(8),intent(out):: mean,var
  integer :: i

  if(size(wgts)/=size(values)) stop 'wgts and values not same size'
  mean = 0d0
  do i=1,size(wgts)
    mean = mean + wgts(i) *  values(i)
  enddo

  var  = 0d0
  do i=1,size(wgts)
    var  = var  + wgts(i) * (values(i)-mean)**2
  enddo
  end subroutine



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



  subroutine mean_var_and_SEM_s( values,numofvalues,mean,var,SEM )
  !This subroutine calculates mean, variance, and standard error of the mean
  !for a sample from a population not the whole population
  real(8) :: values(:),mean,var,SEM
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

  SEM = sqrt(var/numofvalues)

  end subroutine mean_var_and_SEM_s





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



  recursive function HermiteProbpoly(n,x, oin,val_in,val__in) result(val)
  ! This function samples from the probabalists' version of the Hermite polynomials
  ! of order 'n' at location 'x'.  A standard Gaussian is assumed as the basis.
  ! 'o', 'val_', and 'val__' are internal recursion variables that need not be passed when
  ! calling this function.
  integer :: n,   o
  real(8) :: x,val,   val_,val__
  integer,optional :: oin
  real(8),optional :: val_in, val__in

  if(present(oin)) then
    o = oin
  else
    o = 0
  endif
  if(present(val_in))  val_  = val_in
  if(present(val__in)) val__ = val__in
  if(o==0) then
    val  = 1d0
    val_ = 0d0
    val__= 0d0
  elseif(o==1) then
    val  = x
    val__= 0d0
  else
    val = x*val_ - (float(o)-1d0)*val__
  endif

  if(o<n) val = HermiteProbpoly(n,x,o+1,val,val_)
  end function HermiteProbpoly




  recursive function HermitePhyspoly(n,x, oin,val_in,val__in) result(val)
  ! This function samples from the physicists' version of the Hermite polynomials
  ! of order 'n' at location 'x'.  exp(-x^2) is assumed as the basis.
  ! 'o', 'val_', and 'val__' are internal recursion variables that need not be passed when
  ! calling this function.
  integer :: n,   o
  real(8) :: x,val,   val_,val__
  integer,optional :: oin
  real(8),optional :: val_in, val__in

  if(present(oin)) then
    o = oin
  else
    o = 0
  endif
  if(present(val_in))  val_  = val_in
  if(present(val__in)) val__ = val__in
  if(o==0) then
    val  = 1d0
    val_ = 0d0
    val__= 0d0
  elseif(o==1) then
    val  = 2d0 * x
    val__= 0d0
  else
    val = 2d0*x*val_ - 2d0*(float(o)-1d0)*val__
  endif

  if(o<n) val = HermitePhyspoly(n,x,o+1,val,val_)
  end function HermitePhyspoly



  function factorial(n)
  !This function takes the factorial of positive integer 'n'.
  integer :: n,i,factorial
  factorial = 1
  do i=1,n
    factorial = factorial * i
  enddo
  end function factorial


  function nCr(n,r)
  !This function find the binomial coefficient of 'r' and 'k' where n=r+k
  integer :: n,r,i,nCr,numprod,denprod
  numprod = 1
  do i=max(n-r,r)+1,n
    numprod = numprod * i
  enddo
  denprod = 1
  do i=1,min(n-r,r)
    denprod = denprod * i
  enddo
  nCr = numprod/denprod
  end function nCr


  function arithmaticsum(nf,nl,N)
  !This function finds the sum of an arithmatic series of 'N' numbers whose
  !first number is nf and last number is nl.
  !Function assumes all numbers are integers.
  integer :: nf,nl,N,arithmaticsum
  arithmaticsum = ( N * ( nf+nl ) ) / 2 
  end function arithmaticsum


  function geometricsum(base,mult,numterms)
  !This function finds the sum of a geometric series of 'numterms' where the
  !first number is 'base' and each number grows by a factor of 'mult'.
  !Function assumes all numbers (including gap) are integers.
  integer :: base, mult, numterms, i, geometricsum
  geometricsum = 0
  do i=1,numterms
    geometricsum = geometricsum + base*mult**(i-1)
  enddo
  end function geometricsum



  function Heavi(arg)
  !The Heaviside function, if arg less than (or equal to) zero, return 0, else 1
  real(8) :: arg, Heavi
  Heavi = merge(-0.0d0,1.0d0,arg<0.0d0)
  end function Heavi


  function OneGaussrandnum(rand1,rand2)
  !Uses Box-Muller transformation to transform two uniformly-distributed
  !random numbers to be Guassian-distributed random numbers.
  !Returns one of these.  Similar subroutine returns both numbers.
  real(8), intent(in)  :: rand1, rand2
  real(8) :: OneGaussrandnum
  real(8) :: pi  = 3.14159265358979323846d0

  OneGaussrandnum = sqrt(-2.0d0*log(rand1)) * cos(2.0d0*pi*rand2)

  end function OneGaussrandnum


  subroutine TwoGaussrandnums(rand1,rand2,randGauss)
  !Uses Box-Muller transformation to transform two uniformly-distributed
  !random numbers to be Guassian-distributed random numbers.
  !Returns both numbers.  Similar function returns one of two numbers.
  real(8), intent(in)  :: rand1, rand2
  real(8), intent(out) :: randGauss(2)
  real(8) :: pi  = 3.14159265358979323846d0

  randGauss(1) = sqrt(-2.0d0*log(rand1)) * cos(2.0d0*pi*rand2)
  randGauss(2) = sqrt(-2.0d0*log(rand1)) * sin(2.0d0*pi*rand2)

  end subroutine TwoGaussrandnums


  
  function exponentialfit(s,a,lamcsig) result(lamcw)
  !This function prints s, a, and lamcsig (log-normal correlation length) to a text file,
  !executes a python script which fits and exponential curve fit to this log-normal
  !covariance data, and reads in the correlation length of the exponential fit
  real(8) :: s,a,lamcsig,lamcw

  open(unit = 101, file = "auxiliary/expfit/s.a.lamcsig.txt")
  write(101,*) s
  write(101,*) a
  write(101,*) lamcsig
  close(101)

  call system("./auxiliary/expfit/fexpfit.py")
  open(unit = 101, file = "auxiliary/expfit/lamcw.txt")
  read(101, *) lamcw
  close(101)

  end function exponentialfit
 


  subroutine create_cubature(Qs,wgts,nodes)
  !This subroutine accesses a python script which calculates (and arranges in decreasing
  !order of weights) cubature weights and abscissas.
  integer, allocatable :: Qs(:)      !inp
  real(8), allocatable :: wgts(:)    !out
  real(8), allocatable :: nodes(:,:) !out

  integer :: i, q
  real(8), allocatable :: abms(:,:)
  integer :: numpts, numdims
  character(2) :: cubetype = 'GH' !setting to GH now (python code works with GL also)
  print *,"size(nodes(:,1)),size(nodes(1,:)):",size(nodes(:,1)),size(nodes(1,:))
  numdims = size(Qs)
  allocate(abms(numdims,2))

  !hardcode for basic Gaussian, can change later (say offset Gaussian--SCM on Markov-Based)
  do q=1,numdims
    abms(q,1) = 0d0
    abms(q,2) = 1d0
  enddo

  !create input for cubature generation
  open(unit = 101, file = "auxiliary/cubature/Gen_Ord_Cubature.inp.txt")
  write(101,*) cubetype
  write(101,*) numdims
  do i=1,numdims
    write(101,*) abms(i,1),"-",abms(i,2)
  enddo
  1020 format(i5," - ")
  1021 format(i5)
  do i=1,numdims
    if(i<numdims) then
      write(101,1020,advance="no") Qs(i)
    else
      write(101,1021) Qs(i)
    endif
  enddo
  close(101)

  !generate cubature
  call system("./auxiliary/cubature/gen_ordered_cubature.py")

  !read it
  numpts = product(Qs)
  wgts = 0d0
  nodes= 0d0
  open(unit = 101, file = "auxiliary/cubature/Gen_Ord_Cubature.w.x.txt")
  do i=1,numpts
    read(101, *) wgts(i)
  enddo
  print *,"numdims:",numdims
  do i=1,numpts
    do q=1,numdims
      read(101,*) nodes(i,q)
    enddo
  enddo
  close(101)

  end subroutine create_cubature




  subroutine numerical_eigmodesolve(chGausstype,s,procave,procvar,lamc,maxnumeigs,numNystrom,eigs,eigvecs)
  !This subroutine prints values to, executes, and reads values from a python script which solves
  !eigenvalues and eigenvectors for KL expansion.
  character(*) :: chGausstype
  real(8) :: s, procave, procvar, lamc
  integer :: maxnumeigs, numNystrom
  real(8) :: eigs(:),eigvecs(:,:)
  character(3) :: dumchar

  integer :: i,ix

  open(unit = 102, file = "auxiliary/Nystrom/Nystrominp.txt")
  if(chGausstype=='Gaus') then
    write(102,*) 1
  elseif(chGausstype=='LogN') then
    write(102,*) 2
  endif
  write(102,*) s
  write(102,*) procave
  write(102,*) procvar
  write(102,*) lamc
  write(102,*) maxnumeigs
  write(102,*) numNystrom
  close(102)

  call system("./auxiliary/Nystrom/Nystromsolve.py")

  open(unit = 102, file = "auxiliary/Nystrom/Nystromout.txt")
  read(102,*) dumchar
  do i=1,maxnumeigs
    read(102,*) eigs(i)
  enddo
  do i=1,maxnumeigs
    read(102,*) dumchar
    do ix=1,numNystrom
      read(102,*) eigvecs(i,ix)
    enddo
  enddo
  close(102)
  end subroutine numerical_eigmodesolve






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

end module utilities





module mpiaccess 
  implicit none 
  integer :: jobid

#ifdef USE_MPI 
  real(8), dimension(3) :: BATCHINFO 
  integer :: njobs 
  integer :: comm 
  integer, allocatable, dimension(:) :: dutychart
 
  contains 
  subroutine initialize_mpi 
    use mpi 
    implicit none 
    integer :: ierr
    call MPI_INIT(ierr) 
    call MPI_Comm_size(MPI_COMM_WORLD, njobs, ierr)  !determine the number of processes
    call MPI_Comm_rank(MPI_COMM_WORLD, jobid, ierr)  !determine rank of this process
  end subroutine initialize_mpi 

  subroutine assigndutychart_mpi(numpts)
    integer :: i,remainder,numpts
    allocate(dutychart(0:njobs))
    dutychart = numpts/njobs
    dutychart(0) = 0
    remainder = mod(numpts,njobs)
    do i=1,njobs
      if(remainder>0) then 
        dutychart(i) = dutychart(i) + 1 
        remainder = remainder - 1 
      endif 
      if(i>1) dutychart(i) = dutychart(i) + dutychart(i-1)
    enddo 
    if(jobid==njobs-1) print *,"dutychart:",dutychart
  end subroutine assigndutychart_mpi 
#endif 
end module mpiaccess

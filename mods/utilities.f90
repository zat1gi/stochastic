module utilities
  implicit none
CONTAINS
  ! print statements in this module use # 1000-1100

  subroutine calc_time_p( t1,t2,runtime )
  !Calculates runtime in minutes, does print
  real(8) :: t1,t2,runtime
  1000 format("  total runtime:",f10.2," min")

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




!  subroutine store_in_bins( binSmallBound,binLargeBound,binNumof,&
!             binCounts,binBounds,vararray,numVar )
  !Bins values between binSmallBound and binLargeBound with binNumof bins
  !vararray is what you are binning, with numVar as the length of vararray
  !binBounds is bounds with binNumof+1 size, counts in binCounts with binNumof as size
!  integer, intent(in)  :: binNumof,numVar
!  integer, allocatable, intent(out) :: binCounts(:)
!  real(8), intent(in)  :: binSmallBound,binLargeBound,vararray(:)
!  real(8), allocatable, intent(out) :: binBounds(:)

!  integer :: binToStore,j,k
!  real(8) :: binSize

!  allocate(binCounts(binNumof))
!  allocate(binBounds(binNumof+1))

  !set up bins
!  binSize = ( binLargeBound - binSmallBound ) / ( binNumof - 1 )
!  binBounds(1)=binSmallBound
!  do k=2,binNumof+1
!    binBounds(k) = binBounds(k-1) + binSize
!  enddo

  !store in bins
!  do j=1,numVar
!    binToStore = ceiling(( vararray(j) - binSmallBound ) / binSize )
!    binCounts(binToStore) = binCounts(binToStore) + 1
!  enddo

!  deallocate(binCounts)
!  deallocate(binBounds)

!  end subroutine store_in_bins
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
  call system("ps2eps genericplot.ps")  !convert to eps

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




end module utilities

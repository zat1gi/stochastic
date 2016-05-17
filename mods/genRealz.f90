module genRealz
  use utilities
  implicit none

CONTAINS
  ! print statements in this module use # 200-299


  subroutine genbinaryReal( j )
  !creates a binary material realization, plots if specified,
  !and collects tallies for realization stats
  use rngvars, only: rngappnum, rngstride, setrngappnum, mts
  use genRealzvars, only: sig, lam, slen, largesti, numPath, pltgenrealznumof, &
                          nummatSegs, P, matFirstTally, sumPath, sqrPath, &
                          pltgenrealz, matType, matLength, pltgenrealzwhich, totLength
  use mt_stream, only: init, genrand_double3

  integer :: j

  integer, parameter :: numArrSz = 5000 !temp var, don't know how long to make arrays yet
  integer :: i,matType_temp(numArrSz)
  real(8) :: matLength_temp(numArrSz)

  if(allocated(matType)) deallocate(matType)
  if(allocated(matLength)) deallocate(matLength)

  !allows reproducable realizations for same rngseed
  call setrngappnum('genRealz')
  call init( mts , int(rngappnum*rngstride+j,4) )

  matLength_temp=0d0
  matType_temp=0

  if( genrand_double3(mts) < P(1) ) then  !choose which material first and tally
    matType_temp(1) = 1
    matFirstTally(1) = matFirstTally(1)+1
  else
    matType_temp(1) = 2
    matFirstTally(2) = matFirstTally(2)+1
  endif

  matLength_temp(1)=0d0
  i=2
  do while ( matLength_temp(i-1)<slen )

    !decide total length at next segment
    matLength_temp(i)=matLength_temp(i-1)+lam(matType_temp(i-1))*log(1/(1-genrand_double3(mts))) 

    if(matType_temp(i-1)==1) matType_temp(i)=2    !change mat for next seg
    if(matType_temp(i-1)==2) matType_temp(i)=1

    if(matLength_temp(i)>slen) then !truncate if necessary
      matLength_temp(i)=slen
      nummatSegs = i-1
    endif

    numPath(matType_temp(i-1))=numPath(matType_temp(i-1))+1  !tallies for stats

    sumPath(matType_temp(i-1))=sumPath(matType_temp(i-1))+matLength_temp(i)-matLength_temp(i-1)
    sqrPath(matType_temp(i-1))=sqrPath(matType_temp(i-1))+(matLength_temp(i)-matLength_temp(i-1))**2

    if(i>largesti) largesti=1 !track largest i


    if( pltgenrealz(1) .ne. 'noplot' ) then  !print to plot selected realizations
    203 format("          ",A7,"            ",A7)
    205 format(f16.8,"     ",f16.8)

    if( j==pltgenrealzwhich(1) ) then
      if(i==2) open(unit=20, file="genRealzplot1.txt")
      if(i==2) write(20,203) pltgenrealz(3),pltgenrealz(4)
               write(20,205) matLength_temp(i-1),sig(matType_temp(i-1))
               write(20,205) matLength_temp(i),  sig(matType_temp(i-1))
    endif
    if( pltgenrealznumof>=2 .and. j==pltgenrealzwhich(2) ) then
      if(i==2) open(unit=21, file="genRealzplot2.txt")
      if(i==2) write(21,203) pltgenrealz(3),pltgenrealz(4)
               write(21,205) matLength_temp(i-1),sig(matType_temp(i-1))
               write(21,205) matLength_temp(i),  sig(matType_temp(i-1))
    endif
    if( pltgenrealznumof>=3 .and. j==pltgenrealzwhich(3) ) then
      if(i==2) open(unit=22, file="genRealzplot3.txt")
      if(i==2) write(22,203) pltgenrealz(3),pltgenrealz(4)
               write(22,205) matLength_temp(i-1),sig(matType_temp(i-1))
               write(22,205) matLength_temp(i),  sig(matType_temp(i-1))
    endif
    endif

   i=i+1
  enddo

  allocate(matType(nummatSegs))
  allocate(matLength(nummatSegs+1))
  matType   = 0
  matLength = 0.0d0

  do i=1,nummatSegs !translate to allocated arrays
    matType(i)   = matType_temp(i)
    matLength(i) = matLength_temp(i)
  enddo
  matLength(nummatSegs+1) = matLength_temp(nummatSegs+1)

  do i=2,nummatSegs+1 !collect total length data for Actual Co calculations
    totLength(matType(i-1))=totLength(matType(i-1))+matLength(i)-matLength(i-1)
  enddo

  close(unit=20) !close three files printing xs data to
  close(unit=21)
  close(unit=22)

  end subroutine genbinaryReal



  subroutine genLPReal
  !chooses material to come next for Levermore-Pomraning transport
  use genRealzvars, only: P, matType
  use rngvars, only: mts
  use mt_stream, only: genrand_double3

  if(.not.allocated(matType)) allocate(matType(1))
  matType(1) = merge(1,2,genrand_double3(mts)<P(1))

  end subroutine genLPReal




  subroutine genReal_stats
  !stats and plotting for binary material construction
  use genRealzvars, only: sig, lam, numRealz, numPath, pltgenrealznumof, P, &
                          perFirstTally, devFirstTally, matFirstTally, sumPath, &
                          sqrPath, pltgenrealz

  real(8) :: upperbound,lowerbound,avePath(2),devPath(2)

  210 format("Start mat1: ",f8.4,"+-",f8.4,"   mat2: ",f8.4,"+-",f8.4)
  211 format("      P(1): ",f8.4,"             P(2): ",f8.4 )
  212 format("Path  mat1: ",f8.4,"+-",f8.4,"   mat2: ",f8.4,"+-",f8.4)
  213 format("    lam(1): ",f8.4,"           lam(2): ",f8.4 )


  ! ave num of realz that start with certain mat
  perFirstTally(1)=matFirstTally(1)/numRealz
  devFirstTally(1)=sqrt((perFirstTally(1)-perFirstTally(1)**2)/(numRealz-1))
  perFirstTally(2)=matFirstTally(2)/numRealz
  devFirstTally(2)=sqrt((perFirstTally(2)-perFirstTally(2)**2)/(numRealz-1))

  write(*,*)
  write(*,210) perFirstTally(1),devFirstTally(1),perFirstTally(2),devFirstTally(2)
  write(*,211) P(1), P(2)

  ! ave pathlengths yeilded
  avePath(1)=sumPath(1)/numPath(1)
  devPath(1)=sqrt((sqrPath(1)/numPath(1)-avePath(1)**2)/(numPath(1)-1))
  avePath(2)=sumPath(2)/numPath(2)
  devPath(2)=sqrt((sqrPath(2)/numPath(2)-avePath(2)**2)/(numPath(2)-1))

  write(*,*)
  write(*,212) avePath(1),devPath(1),avePath(2),devPath(2)
  write(*,213) lam(1),lam(2)


  !plot if selected
  if(pltgenrealz(1) .NE. 'noplot') then

  open(unit=23, file="genRealzgnuprefix.txt") !create prefix to gnu (for window size)
  214 format("set autoscale")
  215 format("set yrange[",f10.5,":",f10.5,"]")
  upperbound = maxval(sig) + abs( sig(1)-sig(2) ) / 5.00d0
  lowerbound = minval(sig) - abs( sig(1)-sig(2) ) / 5.00d0
  write(23,214)
  write(23,215) lowerbound,upperbound
  close(unit=23)

  if(pltgenrealznumof==1) then  !attach prefix
    call system("cat genRealzgnuprefix.txt plots/genRealz/genRealz.1.gnu > genRealztemp.gnu")
  elseif(pltgenrealznumof==2) then
    call system("cat genRealzgnuprefix.txt plots/genRealz/genRealz.2.gnu > genRealztemp.gnu")
  else
    call system("cat genRealzgnuprefix.txt plots/genRealz/genRealz.3.gnu > genRealztemp.gnu")
  endif

  if(pltgenrealz(1)=='preview') then !attach suffix (pause if desired)
    call system("cat genRealztemp.gnu plots/genRealz/pause.txt > genRealz.gnu")
    call system("rm genRealztemp.gnu")
  else
    call system("mv genRealztemp.gnu genRealz.gnu")
  endif

  call system("gnuplot plots/genRealz/genRealz.gnu") !plot

  call system("ps2pdf genRealzplot.ps") !make pdf
  call system("mv genRealz.gnu genRealzgnuprefix.txt plots/genRealz") !archive
  call system("mv genRealzplot* plots/genRealz")
  endif

  end subroutine genReal_stats





  subroutine reset_genRealtals
  !This subroutine is used when realizations are created more than one time.
  !It resets Markov stats for the next round of tallies
  use genRealzvars, only: numPath, sumPath, sqrPath, largesti, totLength
  numPath   = 0  !setup Markov material tallies
  sumPath   = 0d0
  sqrPath   = 0d0
  largesti  = 0
  totLength = 0d0
  end subroutine reset_genRealtals



end module genRealz

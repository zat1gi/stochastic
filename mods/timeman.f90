module timeman
  implicit none

CONTAINS
  ! print statements in this module use 1100-1200


!! time tracking funcs and subs
  subroutine radtrans_time( j )
  use timevars, only: time
  use genRealzvars, only: numRealz
  use KLvars, only: KLres
  use MCvars, only: trannprt, radMC, radWood
  integer :: j

  real(8) :: timedone,esttime
  character(19) :: type

  if( mod(j,trannprt)==0 ) then
    if(radMC=='yes' .AND. KLres=='yes' .AND. radWood=='yes') type='radMC,KLres,radWood'
    if(radMC=='yes' .AND. KLres=='yes' .AND. radWood=='no')  type='radMC,KLres        '
    if(radMC=='yes' .AND. KLres=='no' .AND. radWood=='yes')  type='radMC,radWood      '
    if(radMC=='yes' .AND. KLres=='no' .AND. radWood=='no')   type='radMC              '
    if(radMC=='no' .AND. KLres=='yes' .AND. radWood=='yes')  type='KLres,radWood      '
    if(radMC=='no' .AND. KLres=='yes' .AND. radWood=='no')   type='KLres              '
    if(radMC=='no' .AND. KLres=='no' .AND. radWood=='yes')   type='radWood            '

    timedone=time(1)
    if(radMC=='yes') timedone=timedone+time(2)
    if(KLres=='yes') timedone=timedone+time(5)
    if(radWood=='yes') timedone=timedone+time(3)
    timedone=timedone/60.0d0
    esttime = timedone*numRealz/j
    call time_report( type,timedone,j,numRealz,esttime )
  endif

  end subroutine radtrans_time




  subroutine KLr_time( j )
  use timevars, only: time
  use KLvars, only: KLrnumRealz
  integer :: j

  real(8) :: timedone,esttime
  character(19) :: type='KLrec              '

  timedone=time(6)/60.0d0
  esttime = timedone*KLrnumRealz/j
  call time_report( type,timedone,j,KLrnumRealz,esttime )

  end subroutine KLr_time



  subroutine KLWood_time( j )
  use timevars, only: time
  use KLvars, only: KLrnumRealz
  integer :: j

  real(8) :: timedone,esttime
  character(19) :: type='KLWood              '

  timedone=time(7)/60.0d0
  esttime = timedone*KLrnumRealz/j
  call time_report( type,timedone,j,KLrnumRealz,esttime )

  end subroutine KLWood_time




  subroutine time_report( type,timedone,j,numj,esttime )
  use timevars, only: t1
  integer :: j,numj
  real(8) :: timedone,esttime
  character(19) :: type

  real(8) :: tottime,t2

  call cpu_time(t2)
  tottime = (t2-t1)/60.0d0

  1101 format(A19,i10," of",i10," realz       tot:",f7.2," time/est:",f7.2,"/",f7.2," min")
  write(*,1101) type,j,numj,tottime,timedone,esttime

  end subroutine time_report




  subroutine timereport
  use timevars, only: t1, runtime
  use KLvars, only: KLres, KLrec, KLnoise
  use MCvars, only: radMC, radWood, KLWood
  use timevars, only: time, ntime, runtime
  use utilities, only: calc_time

  integer :: i
  real(8) :: t2,othertime,otherpercent
  real(8),allocatable :: pertime(:)

  open(unit=100,file="timereport.out")

  call calc_time( t1,t2,runtime )
  119 format("  total runtime:",f10.2," min")
  write(100,119) runtime

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
                     write(100,110) time(1),pertime(1)
  if(radMC=='yes')   write(100,111) time(2),pertime(2)
  if(radWood=='yes') write(100,112) time(3),pertime(3)
  if(KLnoise=='yes') write(100,113) time(4),pertime(4)
  if(KLres=='yes')   write(100,114) time(5),pertime(5)
  if(KLrec=='yes')   write(100,115) time(6),pertime(6)
  if(KLWood=='yes')  write(100,116) time(7),pertime(7)
                     write(100,117) othertime,otherpercent
!  write(100,118) runtime
  close(unit=100)
  call system("mv timereport.out texts")

  end subroutine timereport





end module timeman


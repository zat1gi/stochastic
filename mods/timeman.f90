module timeman
  implicit none

CONTAINS
  ! print statements in this module use 1100-1200


!! time tracking funcs and subs
  subroutine radtrans_time( time,ntime,radMC,KLres,radWood,j,trannprt,t1 )
  use genRealzvars, only: numRealz
  integer :: ntime,j,trannprt
  real(8) :: time(:),t1
  character(3) :: radMC,KLres,radWood

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
    call time_report( type,timedone,j,numRealz,esttime,t1 )
  endif

  end subroutine radtrans_time




  subroutine KLr_time( time,ntime,j,t1 )
  use KLvars, only: KLrnumRealz
  integer :: ntime,j
  real(8) :: time(:),t1

  real(8) :: timedone,esttime
  character(19) :: type='KLrec              '

  timedone=time(6)/60.0d0
  esttime = timedone*KLrnumRealz/j
  call time_report( type,timedone,j,KLrnumRealz,esttime,t1 )

  end subroutine KLr_time



  subroutine KLWood_time( time,ntime,j,t1 )
  use KLvars, only: KLrnumRealz
  integer :: ntime,j
  real(8) :: time(:),t1

  real(8) :: timedone,esttime
  character(19) :: type='KLWood              '

  timedone=time(7)/60.0d0
  esttime = timedone*KLrnumRealz/j
  call time_report( type,timedone,j,KLrnumRealz,esttime,t1 )

  end subroutine KLWood_time




  subroutine time_report( type,timedone,j,numj,esttime,t1 )
  integer :: j,numj
  real(8) :: timedone,esttime,t1
  character(19) :: type

  real(8) :: tottime,t2

  call cpu_time(t2)
  tottime = (t2-t1)/60.0d0

  1101 format(A19,i10," of",i10," realz       tot:",f7.2," time/est:",f7.2,"/",f7.2," min")
  write(*,1101) type,j,numj,tottime,timedone,esttime

  end subroutine time_report


end module timeman


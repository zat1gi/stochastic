module timeman
  implicit none

CONTAINS
  ! print statements in this module use 1100-1200


!! time tracking funcs and subs


  subroutine radtrans_timeupdate( j,icase,tt1 )
  !Print time updates for MCtrans methods.  Also time tracking is performed in part on the 
  !fly in this subroutine.
  use timevars, only: time, totparts, cumparts
  use genRealzvars, only: numRealz
  use MCvars, only: numParts, MCcaseson, MCcases, LPamnumParts, numPosMCmeths
  integer :: j,icase
  real(8) :: tt1

  integer :: ticase
  real(8) :: ttime,localper,timeeta,wgtavetime,tt2
  real(8) :: local_time,finished_time,local_time_left,non_local_time_left
  real(8), allocatable :: avetime(:) !average time per MC method


  !log time (in future subtract out any other contributions above)
  call cpu_time(tt2)
  select case (MCcases(icase))
    case ("radMC")
      time(2)         = time(2)   + (tt2-tt1)
      localper        = real(j,8) / numRealz*100 !percentage of local method done
      cumparts(icase) = j         * numParts     !update cumulative particles
    case ("radWood")
      time(3)         = time(3)   + (tt2-tt1)
      localper        = real(j,8) / numRealz*100
      cumparts(icase) = j         * numParts
    case ("KLWood")
      time(7)         = time(7)   + (tt2-tt1)
      localper        = real(j,8) / numRealz*100
      cumparts(icase) = j         * numParts
    case ("LPMC")
      time(8)         = time(8)   + (tt2-tt1)
      localper        = 100.0d0
      cumparts(icase) = LPamnumParts
    case ("atmixMC")
      time(9)         = time(9)   + (tt2-tt1)
      localper        = 100.0d0
      cumparts(icase) = LPamnumParts
  end select
  tt1 = tt2                         !reset tt1 to tt2

  !get time estimates
  if(.not.allocated(avetime)) allocate(avetime(numPosMCmeths))
  avetime    = 0.0d0
  wgtavetime = 0.0d0
  do ticase=1,icase
    if(MCcaseson(ticase)==1) then
      select case (MCcases(ticase)) !load to ttime the time of each method
        case ("radMC")
          avetime(ticase) = time(2) / cumparts(ticase) !ave time per particle for method
          wgtavetime      = wgtavetime + avetime(ticase) * numparts * numRealz
        case ("radWood")
          avetime(ticase) = time(3) / cumparts(ticase)
          wgtavetime      = wgtavetime + avetime(ticase) * numparts * numRealz
        case ("KLWood")
          avetime(ticase) = time(7) / cumparts(ticase)
          wgtavetime      = wgtavetime + avetime(ticase) * numparts * numRealz
        case ("LPMC")
          avetime(ticase) = time(8) / cumparts(ticase)
          wgtavetime      = wgtavetime + avetime(ticase) * LPamnumparts
        case ("atmixMC")
          avetime(ticase) = time(9) / cumparts(ticase)
          wgtavetime      = wgtavetime + avetime(ticase) * LPamnumparts
      end select
    endif
  enddo
  wgtavetime = wgtavetime / sum(totparts) !weighted average time per particle estimate

  local_time          = cumparts(icase)                                                  *avetime(icase)

  finished_time       = sum(time)
  local_time_left     = (totparts(icase)-cumparts(icase))                                *avetime(icase)
  non_local_time_left = (sum(totparts)-sum(cumparts) - (totparts(icase)-cumparts(icase)))*wgtavetime
!print *,"local_time_left: ",local_time_left," non_local_time_left:",non_local_time_left
  timeeta             = finished_time + local_time_left + non_local_time_left


  1100 format(A9,"   ",f7.2," min,",f6.1,"% of method ",i2," of ",i2,"      time/est:",f7.2,"/",f7.2," min")
  write(*,1100) MCcases(icase),local_time/60.0d0,localper,sum(MCcaseson(1:icase)),sum(MCcaseson),&
                finished_time/60.0d0,timeeta/60.0d0
  if(localper==100) print *,

  end subroutine radtrans_timeupdate





  subroutine KL_timeupdate( j,tt1,flKLtype )
  !Print time updates for KLcol and KLrec methods.  Also time tracking is performed in part on the 
  !fly in this subroutine.
  use timevars, only: time
  use genRealzvars, only: numRealz
  use KLvars, only: KLrnumRealz
  integer :: j
  real(8) :: tt1,localper,local_time,timeeta
  character(5) :: flKLtype

  real(8) :: tt2


  call cpu_time(tt2)     !increment time
  select case (flKLtype)
    case ("KLcol")
      time(5) = time(5) + (tt2-tt1)
    case ("KLrec")
      time(6) = time(6) + (tt2-tt1)
  end select      
  tt1 = tt2


  1110 format(A9,"   ",f7.2," min,",f6.1,"% of est",f7.2," min for meth")
  select case (flKLtype) !print update
    case ("KLcol")
      localper   = real(j,8)/numRealz*100
      local_time = time(5)
      timeeta    = local_time / (localper/100)
      write(*,1110) flKLtype,local_time/60.0d0,localper,timeeta/60.0d0
    case ("KLrec")
      localper   = real(j,8)/KLrnumRealz*100
      local_time = time(6)
      timeeta    = local_time / (localper/100)
      write(*,1110) flKLtype,local_time/60.0d0,localper,timeeta/60.0d0
  end select      


  if(flKLtype=='KLcol' .and. j==numRealz) time(5) = time(5)

  end subroutine KL_timeupdate





  subroutine timereport
  !This is the timereport that comes at the end of the whole code running, the final time report.
  use timevars, only: t1, runtime
  use KLvars, only: KLres, KLrec, KLnoise
  use MCvars, only: radMC, radWood, KLWood, LPMC, atmixMC
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
  othertime = runtime + time(1) / 60 !'time(1)' is genRealz.  It is incorporated within other times already.
                                !This addition compensates for otherwise double counting.
  do i=1,ntime
    time(i)    = time(i) / 60
    othertime  = othertime - time(i)
    pertime(i) = time(i) / runtime * 100
  enddo
  otherpercent = othertime / runtime * 100

  110 format("    genRealz   :",f10.2," min     per:",f6.2," % --- genRealz in other meths,")
  121 format("    ----------------------------------------------- not part of sum")
  111 format("    radMC      :",f10.2," min     per:",f6.2," %")
  112 format("    radWood    :",f10.2," min     per:",f6.2," %")
  113 format("    KLnoise    :",f10.2," min     per:",f6.2," %")
  114 format("    KLcol      :",f10.2," min     per:",f6.2," %")
  115 format("    KLrec      :",f10.2," min     per:",f6.2," %")
  116 format("    KLWood     :",f10.2," min     per:",f6.2," %")
  117 format("    other      :",f10.2," min     per:",f6.2," %")
  118 format("    LPMC       :",f10.2," min     per:",f6.2," %")
  120 format("    atmixMC    :",f10.2," min     per:",f6.2," %")

                     write(100,110) time(1),pertime(1)
                     write(100,121)
  if(radMC=='yes')   write(100,111) time(2),pertime(2)
  if(radWood=='yes') write(100,112) time(3),pertime(3)
  if(KLWood=='yes')  write(100,116) time(7),pertime(7)
  if(LPMC=='yes')    write(100,118) time(8),pertime(8)
  if(atmixMC=='yes') write(100,120) time(9),pertime(9)
  if(KLnoise=='yes') write(100,113) time(4),pertime(4)
  if(KLres=='yes')   write(100,114) time(5),pertime(5)
  if(KLrec=='yes')   write(100,115) time(6),pertime(6)
                     write(100,117) othertime,otherpercent
!  write(100,118) runtime
  close(unit=100)
  call system("mv timereport.out texts")

  end subroutine timereport










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



end module timeman


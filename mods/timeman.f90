module timeman
  implicit none

CONTAINS
  ! print statements in this module use 1100-1200


!! time tracking funcs and subs


  subroutine radtrans_timeupdate( j,tt1 )
  !Print time updates for MCtrans methods.  Also time tracking is performed in part on the 
  !fly in this subroutine.
  use timevars, only: time, totparts, cumparts
  use genRealzvars, only: numRealz
  use MCvars, only: numParts, chTrantype, LPamnumParts
  integer :: j
  real(8) :: tt1

  real(8) :: ttime,localper,timeeta,wgtavetime,tt2
  real(8) :: local_time,finished_time,local_time_left,non_local_time_left
  real(8), allocatable :: avetime !average time per MC method


  !log time (in future subtract out any other contributions above)
  call cpu_time(tt2)
  select case (chTrantype)
    case ("radMC")
      time(2)         = time(2)   + (tt2-tt1)
      localper        = real(j,8) / numRealz*100 !percentage of local method done
      cumparts = j         * numParts     !update cumulative particles
    case ("radWood")
      time(3)         = time(3)   + (tt2-tt1)
      localper        = real(j,8) / numRealz*100
      cumparts = j         * numParts
    case ("KLWood")
      time(7)         = time(7)   + (tt2-tt1)
      localper        = real(j,8) / numRealz*100
      cumparts = j         * numParts
    case ("LPMC")
      time(8)         = time(8)   + (tt2-tt1)
      localper        = 100.0d0
      cumparts = LPamnumParts
    case ("atmixMC")
      time(9)         = time(9)   + (tt2-tt1)
      localper        = 100.0d0
      cumparts = LPamnumParts
    case ("GaussKL")
      time(10)        = time(10)  + (tt2-tt1)
      localper        = real(j,8) / numRealz*100
      cumparts = j         * numParts
  end select
  tt1 = tt2                         !reset tt1 to tt2

  !get time estimates
  avetime    = 0.0d0
  wgtavetime = 0.0d0
  select case (chTrantype) !load to ttime the time of each method
    case ("radMC")
      avetime = time(2) / cumparts !ave time per particle for method
      wgtavetime      = wgtavetime + avetime * numparts * numRealz
    case ("radWood")
      avetime = time(3) / cumparts
      wgtavetime      = wgtavetime + avetime * numparts * numRealz
    case ("KLWood")
      avetime = time(7) / cumparts
      wgtavetime      = wgtavetime + avetime * numparts * numRealz
    case ("LPMC")
      avetime = time(8) / cumparts
      wgtavetime      = wgtavetime + avetime * LPamnumparts
    case ("atmixMC")
      avetime = time(9) / cumparts
      wgtavetime      = wgtavetime + avetime * LPamnumparts
    case ("GaussKL")
      avetime = time(10)/ cumparts
      wgtavetime      = wgtavetime + avetime * numparts * numRealz
  end select
  wgtavetime = wgtavetime / sum(totparts) !weighted average time per particle estimate

  local_time          = cumparts                                                  *avetime

  finished_time       = sum(time)
  local_time_left     = (totparts-cumparts)                                *avetime
  non_local_time_left = (sum(totparts)-sum(cumparts) - (totparts-cumparts))*wgtavetime
!print *,"local_time_left: ",local_time_left," non_local_time_left:",non_local_time_left
  timeeta             = finished_time + local_time_left + non_local_time_left


  1100 format(A9,"   ",f7.2," min,",f6.1,"% of method   time/est:",f7.2,"/",f7.2," min")
  write(*,1100) chTrantype,local_time/60.0d0,localper,finished_time/60.0d0,timeeta/60.0d0
  if(localper==100) print *,

  end subroutine radtrans_timeupdate





  subroutine KL_timeupdate( j,tt1,flKLtype )
  !Print time updates for KLcol and KLrec methods.  Also time tracking is performed in part on the 
  !fly in this subroutine.
  use timevars, only: time
  use genRealzvars, only: numRealz
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
      localper   = real(j,8)/numRealz*100
      local_time = time(6)
      timeeta    = local_time / (localper/100)
      write(*,1110) flKLtype,local_time/60.0d0,localper,timeeta/60.0d0
  end select      


  if(flKLtype=='KLcol' .and. j==numRealz) time(5) = time(5)

  end subroutine KL_timeupdate





  subroutine timereport
  !This is the timereport that comes at the end of the whole code running, the final time report.
  use timevars, only: t1, runtime, time, ntime, runtime, FOM
  use KLvars, only: KLres, KLrec, KLnoise
  use MCvars, only: stocMC_transmission, stocMC_reflection
  use utilities, only: calc_time

  integer :: i
  real(8) :: t2,othertime,otherpercent
  real(8),allocatable :: pertime(:)

  open(unit=100,file="timereport.out")

  call calc_time( t1,t2,runtime )
  119 format("  total runtime:",f10.2," min")
  write(100,119) runtime

  allocate(pertime(ntime))
  othertime = runtime + time(1) / 60 + time(6) / 60 !'time(1)' is genRealz, 'time(6) is KLrec
                                                    !These are incorporated within other times already.
                                                    !This addition compensates for otherwise double counting.
  do i=1,ntime
    time(i)    = time(i) / 60
    othertime  = othertime - time(i)
    pertime(i) = time(i) / runtime * 100
  enddo
  otherpercent = othertime / runtime * 100

  110 format("    genRealz   :",f7.2," min   per:",f6.2," % --- genRealz in others, not part of sum")
  115 format("    KLrec      :",f7.2," min   per:",f6.2," % ---     Krec in others, not part of sum")
  121 format("    ------------------------------------------------ refl ----- tran --")
  111 format("    radMC      :",f7.2," min   per:",f6.2," %   FOM:  ",es9.3,"  ",es9.3)
  112 format("    radWood    :",f7.2," min   per:",f6.2," %   FOM:  ",es9.3,"  ",es9.3)
  113 format("    KLnoise    :",f7.2," min   per:",f6.2," %")
  114 format("    KLcol      :",f7.2," min   per:",f6.2," %")
  116 format("    KLWood     :",f7.2," min   per:",f6.2," %   FOM:  ",es9.3,"  ",es9.3)
  117 format("    other      :",f7.2," min   per:",f6.2," %")
  118 format("    LPMC       :",f7.2," min   per:",f6.2," %")
  120 format("    atmixMC    :",f7.2," min   per:",f6.2," %")
  122 format("    GaussKL    :",f7.2," min   per:",f6.2," %   FOM:  ",es9.3,"  ",es9.3)


                     write(100,110) time(1), pertime(1)
  if(KLrec=='yes')   write(100,115) time(6), pertime(6)
                     write(100,121)
  if(chTrantype=='radMC') write(100,111) time(2), pertime(2), FOM(1,1),FOM(1,2)
  if(chTrantype=='radWood') write(100,112) time(3), pertime(3), FOM(2,1),FOM(2,2)
  if(chTrantype=='KLWood')  write(100,116) time(7), pertime(7), FOM(3,1),FOM(3,2)
  if(chTrantype=='LPMC')    write(100,118) time(8), pertime(8)
  if(chTrantype=='atmixMC') write(100,120) time(9), pertime(9)
  if(chTrantype=='GaussKL') write(100,122) time(10),pertime(10),FOM(7,1),FOM(7,2)
  if(KLnoise=='yes') write(100,113) time(4), pertime(4)
  if(KLres=='yes')   write(100,114) time(5), pertime(5)
                     write(100,117) othertime,otherpercent
!  write(100,118) runtime
  close(unit=100)
  call system("mv timereport.out texts")

  end subroutine timereport



  subroutine FOM_calculation
  !This subroutine calculates a FOM for all realization-based cases.
  !The FOM = 1/(R^2*t) where R is relative error and t is runtime.
  !Said another way FOM = (ave^2*N)/(stdev^2*t).
  use genRealzvars, only: numRealz
  use MCvars, only: chTrantype, stocMC_reflection, stocMC_transmission
  use timevars, only: FOM, time

  select case (chTrantype)
    case ("radMC")
      FOM(1)=stocMC_reflection(1)**2*real(numRealz,8)/(stocMC_reflection(2)**2*time(2))
      FOM(2)=stocMC_transmission(1)**2*real(numRealz,8)/(stocMC_transmission(2)**2*time(2))
    case ("radWood")
      FOM(1) = stocMC_reflection(1)**2*real(numRealz,8)/(stocMC_reflection(2)**2*time(3))
      FOM(2) = stocMC_transmission(1)**2*real(numRealz,8)/(stocMC_transmission(2)**2*time(3))
    case ("KLWood")
      FOM(1)=stocMC_reflection(1)**2*real(numRealz,8)/(stocMC_reflection(2)**2*time(7))
      FOM(2)=stocMC_transmission(1)**2*real(numRealz,8)/(stocMC_transmission(2)**2*time(7))
    case ("LPMC")
    case ("atmixMC")
    case ("GaussKL")
      FOM(1)=stocMC_reflection(1)**2*real(numRealz,8)/(stocMC_reflection(2)**2*time(10))
      FOM(2)=stocMC_transmission(1)**2*real(numRealz,8)/(stocMC_transmission(2)**2*time(10))
  end select

  end subroutine FOM_calculation



end module timeman


module timeman
  implicit none

CONTAINS
  ! print statements in this module use 1100-1200


!! time tracking funcs and subs

  subroutine initialize_t1
  !sets t1, the beginning of a new repetative operation
  use timevars, only: t1
  call cpu_time(t1)
  end subroutine initialize_t1



  subroutine timeupdate( chtype,ndone,ntotal )
  !uses t1 and new t2 to give update on time for new repetative operation
  use timevars, only: t1
  integer :: ndone,ntotal
  real(8) :: t2
  character(*) :: chtype

  call cpu_time(t2)

  1100 format(A9,"   ",f6.1,"% of ",i7,"   time/est:",f7.2,"/",f7.2," min")
  write(*,1100) chtype,real(ndone,8)/ntotal*100d0,ntotal,(t2-t1)/60.0d0,(t2-t1)*ntotal/ndone/60.0d0
  end subroutine


end module timeman


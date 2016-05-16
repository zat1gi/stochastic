program mt_stream_test
  use mt_stream, only: set_mt19937, new, init, mt_state, genrand_double3
  implicit none
  ! Mersenne Twister testing main

  type(mt_state), save :: mts
  real(8) :: rand
  integer :: iseed1 = 1234
  integer :: iseed2 = 1235

  print *,"This is before everything"
  call set_mt19937
  call new(mts)

  call init(mts,iseed1)
  rand = genrand_double3(mts)
  print *,"rand:",rand
  rand = genrand_double3(mts)
  print *,"rand:",rand

  call init(mts,iseed2)
  rand = genrand_double3(mts)
  print *,"rand:",rand

  call init(mts,iseed1)
  rand = genrand_double3(mts)
  print *,"rand:",rand
  rand = genrand_double3(mts)
  print *,"rand:",rand

  print *,"And this is after"

end program mt_stream_test

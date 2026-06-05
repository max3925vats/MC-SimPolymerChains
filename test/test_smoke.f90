program test_smoke
  use polymc_kinds, only: dp
  implicit none
  if (kind(1.0_dp) /= kind(1.0d0)) stop 1
  print *, "smoke ok"
end program test_smoke

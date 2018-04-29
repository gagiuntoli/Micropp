program main_wrap

  integer dims
  integer :: sizes(3)
  real*8 :: strain(3)
  real*8 :: stress(3)
  real*8 :: ctan(9)
  integer cg_its
  real*8 cg_tol

  dims = 2 
  sizes(1) = 10
  sizes(2) = 10
  sizes(3) = 1
  cg_its = 999
  cg_tol = +1.0e-8

  call micro_construct(dims, sizes, cg_its, cg_tol) 

  strain(1) = 0.005
  strain(2) = 0.0
  strain(3) = 0.0
  call micro_loc_hom_stress(strain, stress) 
  WRITE(*,'(F12.2,F12.2,F12.2,A)') stress(1), stress(2), stress(3)

  call micro_loc_hom_ctan(strain, ctan) 
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(1), ctan(2), ctan(3)
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(4), ctan(5), ctan(6)
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(7), ctan(8), ctan(9)

end program

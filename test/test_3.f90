program main_wrap

  implicit none

  integer dims
  integer :: sizes(3), types(2)
  real*8 :: strain(3)
  real*8 :: stress(3)
  real*8 :: ctan(9)
  integer micro_type
  real*8 :: micro_params(4), params(20)

  dims = 2 
  sizes(1) = 10
  sizes(2) = 10
  sizes(3) = 1

  micro_type = 0; ! 2 materiales matriz y fibra (3D esfera en matriz)
  micro_params(1) = 1.0; ! lx
  micro_params(2) = 1.0; ! ly
  micro_params(3) = 1.0; ! lz
  micro_params(4) = 0.2; ! radio de la esfera
  types(1) = 0;
  types(2) = 0;

  params(0*10 + 1) = 1.0e6;
  params(0*10 + 2) = 0.3;
  params(1*10 + 1) = 1.0e7;
  params(1*10 + 2) = 0.3;

  call micro_construct(dims, sizes, micro_type, micro_params, types, params)

  strain(1) = 0.005
  strain(2) = 0.0
  strain(3) = 0.0
  call micro_loc_hom_stress(69,strain, stress) 
  WRITE(*,'(F12.2,F12.2,F12.2,A)') stress(1), stress(2), stress(3)

  call micro_loc_hom_ctan(69,strain, ctan) 
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(1), ctan(2), ctan(3)
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(4), ctan(5), ctan(6)
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(7), ctan(8), ctan(9)

end program

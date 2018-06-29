!
!  MicroPP : 
!  Finite element library to solve microstructural problems for composite materials.
!
! Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!

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
  micro_params(4) = 0.0; ! I max
  types(1) = 0;
  types(2) = 0;

  params(0*10 + 1) = 1.0e6;
  params(0*10 + 2) = 0.3;
  params(1*10 + 1) = 1.0e7;
  params(1*10 + 2) = 0.3;

  call micropp_construct(dims, sizes, micro_type, micro_params, types, params)

  strain(1) = 0.005
  strain(2) = 0.0
  strain(3) = 0.0
  call micropp_set_macro_strain(69, strain) 
  call micropp_localize_homogenize()
  call micropp_get_macro_stress(69, stress) 

  WRITE(*,*) 'Stress : '
  WRITE(*,'(F12.2,F12.2,F12.2,A)') stress(1), stress(2), stress(3)

  call micropp_get_macro_ctan(69, ctan) 
  WRITE(*,*) 'Ctan : '
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(1), ctan(2), ctan(3)
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(4), ctan(5), ctan(6)
  WRITE(*,'(F12.2,F12.2,F12.2,A)') ctan(7), ctan(8), ctan(9)

end program

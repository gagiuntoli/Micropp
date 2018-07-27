!
! This is a test example for MicroPP: a finite element library
! to solve microstructural problems for composite materials.
!
!  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
!                         Guido Giuntoli <gagiuntoli@gmail.com>
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

program test3d_3

  ! access computing environment
  use ISO_FORTRAN_ENV, only : ERROR_UNIT
  use libmicropp

  implicit none

  integer, parameter :: MAX_MAT_PARAM = 10
  type(micropp3) :: micro
  integer :: argc, t
  character(len=32) :: arg
  integer :: sizes(3), time_steps

  integer, parameter :: micro_type = 1
  real(8), parameter :: d_eps = 0.01
  integer, parameter :: dir = 3;

  real(8), dimension(*) :: eps(6), sig(6)

  real(8) :: micro_params(5)
  type(material_t) :: mat_params(2)

  argc = command_argument_count()

  if (argc < 3) then
     call get_command_argument(0, arg)
     write (ERROR_UNIT,*) "Usage: ./executable nx ny nz [steps]"
     stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) sizes(1)
  call get_command_argument(2, arg)
  read(arg, *) sizes(2)
  call get_command_argument(3, arg)
  read(arg, *) sizes(3)

  if (argc >= 4) then
     call get_command_argument(4, arg)
     read(arg, *) time_steps
  else
     time_steps = 10
  end if

  micro_params = (/ 1.0, 1.0, 1.0, 0.1, 0.0 /)

  call set(mat_params(1), 1.0e6, 0.3, 5.0e4, 5.0e4, 1)
  call set(mat_params(2), 1.0e6, 0.3, 1.0e4, 0.0e-1, 0)

  micro = micropp3(1, sizes, micro_type, micro_params, mat_params)

  eps = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  sig = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  do t = 0, time_steps - 1
     write (*,*) "time step = ", t
     if (t < 30) then
        eps(dir) = eps(dir) + d_eps;
     else if (t < 80) then
        eps(dir) = eps(dir) - d_eps;
     else if (t < 130) then
        eps(dir) = eps(dir) + d_eps;
     else if (t < 250) then
        eps(dir) = eps(dir) - d_eps;
     else
        eps(dir) = eps(dir) + d_eps;
     end if

     call micro%set_macro_strain(0, eps)
     call micro%homogenize()

     call micro%get_macro_stress(0, sig)

     call micro%update_vars()
     call micro%write_info_files();

     write(*,'(A,F12.2)') "eps = ", eps(dir)
     write(*,'(A)', advance="no") 'sig = '
     write(*,'(F12.2,F12.2,F12.2,A)', advance="no") sig(1), sig(2), sig(3)
     write(*,'(F12.2,F12.2,F12.2,A)') sig(4), sig(5), sig(6)
     call micro%output(t, 0);
  end do

  call free(micro)

end program test3d_3

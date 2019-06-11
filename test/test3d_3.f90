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

        type(micropp3) :: micro
        integer :: argc, t
        integer :: i, j
        character(len=32) :: arg
        integer :: sizes(3), time_steps
        logical :: non_linear, converged, subiterated

        integer, parameter :: gp_id = 0
        integer, parameter :: micro_type = 2
        real(8), parameter :: d_eps = 0.01
        integer, parameter :: dir = 3;
        integer :: cost
        integer :: nn
        integer :: nsubiterations
        integer :: mpi_rank
        Character(len = 128) :: filename
        Character(len = 16) :: time_char

        real(8), dimension(*) :: eps(6), sig(6), ctan(36)
        integer, dimension(:) :: coupling(1)

        real(8) :: micro_params(4)
        type(material_base) :: mat_params(3)

        argc = command_argument_count()

        if (argc < 1) then
                call get_command_argument(0, arg)
                write (ERROR_UNIT,*) "Usage: ./executable n [steps(10)]"
                stop 1
        end if

        call get_command_argument(1, arg)
        read(arg, *) nn

        sizes = (/ nn, nn, nn /)

        if (argc >= 2) then
                call get_command_argument(2, arg)
                read(arg, *) time_steps
        else
                time_steps = 10
        end if

        nsubiterations = 1
        mpi_rank = 0
        micro_params = (/ 1.0, 1.0, 1.0, 0.1 /)
        coupling(1) = 1

        call material_set(mat_params(1), 1, 1.0D6, 0.3D0, 5.0D4, 5.0D4, 3.0D3)
        call material_set(mat_params(2), 1, 1.0D6, 0.3D0, 1.0D4, 0.0D1, 3.0D3)
        call material_set(mat_params(3), 1, 1.0D6, 0.3D0, 1.0D4, 0.0D1, 3.0D3)

        call micropp3_new(micro, 1, sizes, micro_type, micro_params, &
                mat_params, coupling, nsubiterations, mpi_rank)
        call micropp3_print_info(micro)

        eps = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
        sig = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
        do t = 0, time_steps - 1
        write (*,'(A,I5)') "time step = ", t
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

        call micropp3_set_strain(micro, gp_id, eps)
        call micropp3_homogenize(micro)

        call micropp3_get_stress(micro, gp_id, sig)
        call micropp3_get_ctan(micro, gp_id, ctan)

        non_linear = micropp3_is_non_linear(micro, gp_id)
        converged = micropp3_has_converged(micro, gp_id)
        subiterated = micropp3_has_subiterated(micro, gp_id)
        cost = micropp3_get_cost(micro, gp_id)

        call micropp3_update_vars(micro)

        write(*,'(A,L)') "Non-linear = ", non_linear
        write(*,'(A,2I5)') "Cost       = ", cost
        write(*,'(A,L)') "Converged  = ", converged
        write(*,'(A,L)') "Subiterated  = ", subiterated
        write(*,'(A,F12.2)') "eps = ", eps(dir)

        write(*,'(A)', advance="no") 'sig = '
        do i = 1, 6
                write(*,'(F12.2)', advance="no") sig(i)
        enddo
        write(*,*)

        write(*,'(A)', advance="no") 'ctan = '
        write(*,*)
        do i = 0, 5
                do j = 1, 6
                        write(*,'(F12.2)', advance="no") ctan(i * 6 + j)
                enddo
                write(*,*)
        enddo

        write(*,*) ""

        !write(filename, "(A,I1)") "micropp_fortran_", t
        call micropp3_output2(micro, 0, 1234, t) 

        end do

        call micropp3_free(micro)

end program test3d_3

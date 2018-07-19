!
! This file is part of the Fortran_C++ program.
! Copyright (c) 2018 Jimmy Aguilar Mena.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!

module libmicropp
  implicit none

  private :: new_micropp3, set_macro_strain, homogenize, &
       get_macro_stress, update_vars, write_info_files

  public :: micropp3, free

  type :: micropp3
     private
     integer(8) :: ptr ! pointer
   contains
     procedure :: print_info, set_macro_strain, homogenize,&
          get_macro_stress, update_vars, write_info_files
  end type micropp3

  interface micropp3
     procedure new_micropp3
  end interface micropp3

  integer(8) :: init3
  external init3

contains
  function new_micropp3(ngp, size, micro_type, micro_params,&
       mat_types, params)
    implicit none
    type(micropp3) :: new_micropp3
    integer, intent(in) :: ngp
    integer, intent(in) :: size(3)
    integer, intent(in) :: micro_type
    real(8), intent(in), dimension (*) :: micro_params
    integer, intent(in), dimension (*) :: mat_types
    real(8), intent(in), dimension (*) :: params

    new_micropp3%ptr = init3(ngp, size, micro_type, micro_params, &
         mat_types, params)

  end function new_micropp3

  subroutine free(this)
    class(micropp3) :: this
    call free3(this%ptr)
  end subroutine free

  subroutine set_macro_strain(this, gp_id, macro_strain)
    class(micropp3), intent(inout) :: this
    integer, intent(in) :: gp_id
    real(8), intent(in), dimension(*) :: macro_strain
    call set_macro_strain3(this%ptr, gp_id, macro_strain)
  end subroutine set_macro_strain

  subroutine homogenize(this)
    class(micropp3) :: this
    call homogenize3(this%ptr)
  end subroutine homogenize

  subroutine get_macro_stress(this, gp_id, macro_stress)
    class(micropp3) :: this
    integer, intent(in) :: gp_id
    real(8), intent(in), dimension(*) :: macro_stress
    call get_macro_stress3(this%ptr, gp_id, macro_stress)
  end subroutine get_macro_stress

  subroutine update_vars(this)
    class(micropp3) :: this
    call update_vars3(this%ptr)
  end subroutine update_vars

  subroutine write_info_files(this)
    class(micropp3) :: this
    call write_info_files3(this%ptr)
  end subroutine write_info_files

  subroutine print_info(this)
    class(micropp3) :: this
    call print_info3(this%ptr)
  end subroutine print_info

end module libmicropp

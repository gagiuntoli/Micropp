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

module libmaterial
  use iso_c_binding
  implicit none

  public :: material_t

  type, bind(C) :: material_t
     real(c_double) :: E, nu, Ka, Sy
     real(c_double) :: k, mu, lambda
     integer(c_int) type
     logical(c_bool) plasticity, damage
  end type material_t

contains

  subroutine set(this, E, nu, Ka, Sy, type)
    implicit none
    type(material_t) :: this
    real(4), intent(in) :: E, nu, Ka, Sy
    integer(4), intent(in) :: type

    call material_set(this, dble(E), dble(nu), & 
         dble(Ka), dble(Sy), type)

  end subroutine set

end module libmaterial

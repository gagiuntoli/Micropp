/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#pragma once


#ifdef _CUDA
#define CUDA_HOSTDEV __host__ __device__

struct cuda_params_t {

	int *elem_type_d;

};

#else
#define CUDA_HOSTDEV
#endif


#define DIM 3
#define NPE 8
#define NVOI 6


CUDA_HOSTDEV
#pragma acc routine seq
void get_elem_nodes(int n[8], const int nx, const int ny,
		    const int ex, const int ey, const int ez = 0);

CUDA_HOSTDEV
#pragma acc routine seq
void get_elem_displ(const double *u, double elem_disp[NPE * DIM],
		    int nx, int ny, int ex, int ey, int ez);

CUDA_HOSTDEV
#pragma acc routine seq
void get_strain(const double *u, int gp, double *strain_gp,
		const double bmat[NPE][NVOI][NPE * DIM], 
		int nx, int ny, int ex, int ey, int ez);

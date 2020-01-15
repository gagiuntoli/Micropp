/*
 *  This source code is part of Micropp: a Finite Element library
 *  to solve composite materials micro-scale problems.
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


#include "micropp.hpp"
#include "common.hpp"
#include "cuda.hpp"


__global__
void device_init_material()
{
	for (int i = 0; i < MAX_MATERIALS; ++i) {
		material_list_d[i] = new material_elastic(1.0, 1.0);
	}
}

__global__
void device_delete_material()
{
	for (int i = 0; i < MAX_MATERIALS; ++i) {
		delete material_list_d;
	}
}

template<>
void micropp<3>::cuda_init()
{
	cudaMalloc((void **)&cuda_params.elem_type_d, nelem * sizeof(int));

	device_init_material<<<1, 1>>>();

	cudaMemcpy(cuda_params.elem_type_d, elem_type, 
		   nelem * sizeof(int), cudaMemcpyHostToDevice);
}

template<>
void micropp<3>::cuda_finalize()
{
	cudaFree(cuda_params.elem_type_d);
}

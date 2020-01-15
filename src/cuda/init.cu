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

__device__ material_t *material_list_d[MAX_MATERIALS];

__global__
void device_init_material()
{
	//for (int i = 0; i < MAX_MATERIALS; ++i) {
	//	material_list_d[i] = new material_elastic(1.0, 1.0);
	//}
}

template<>
void micropp<3>::cuda_init()
{
	cudaMalloc((void **)&cuda_params.elem_type_d, nelem * sizeof(int));
	for (int i = 0; i < MAX_MATERIALS; ++i) {
		cudaMalloc((void **)&cuda_params.material_list, 
				MAX_MATERIALS * sizeof(material_t));
	}

	cudaMemcpy(cuda_params.elem_type_d, elem_type, 
		   nelem * sizeof(int), cudaMemcpyHostToDevice);
	for (int i = 0; i < MAX_MATERIALS; ++i) {
		cudaMemcpy(&cuda_params.material_list[i], &material_list[i], 
			   sizeof(material_t), cudaMemcpyHostToDevice);
	}
}

template<>
void micropp<3>::cuda_finalize()
{
	cudaFree(cuda_params.elem_type_d);
}

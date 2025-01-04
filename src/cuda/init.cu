/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "common.hpp"
#include "cuda.hpp"
#include "micropp.hpp"

__device__ material_t *material_list_d[MAX_MATERIALS];
material_base *material_base_list_d;
// struct cuda_params_t cuda_params;
int *elem_type_d;

__global__ void device_init_material(material_base *material_base_list_d) {
  for (int i = 0; i < MAX_MATERIALS; ++i) {
    material_list_d[i] = material_t::make_material(material_base_list_d[i]);
  }
}

__global__ void device_delete_material() {
  for (int i = 0; i < MAX_MATERIALS; ++i) {
    delete material_list_d;
  }
}

template <>
void micropp<3>::cuda_init(const micropp_params_t &params) {
  cudaMalloc((void **)&elem_type_d, nelem * sizeof(int));

  cudaMalloc((void **)&material_base_list_d, MAX_MATERIALS * sizeof(material_base));
  for (int i = 0; i < MAX_MATERIALS; ++i) {
    cudaMemcpy(&material_base_list_d[i], &params.materials[i], sizeof(material_base), cudaMemcpyHostToDevice);
  }

  device_init_material<<<1, 1>>>(material_base_list_d);

  cudaMemcpy(elem_type_d, elem_type, nelem * sizeof(int), cudaMemcpyHostToDevice);
}

template <>
void micropp<3>::cuda_finalize() {
  // cudaFree(cuda_params.elem_type_d);
}

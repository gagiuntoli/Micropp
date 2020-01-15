
#pragma once


#include "params.hpp"

__device__ material_t *material_list_d[MAX_MATERIALS];

struct cuda_params_t {

	int *elem_type_d;

};

struct cuda_params_t cuda_params;

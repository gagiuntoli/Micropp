
#include "micropp.hpp"
#include "common.hpp"

template<>
void micropp<3>::cuda_init()
{
	cudaMalloc((void **)cuda_params.elem_type_d, nelem * sizeof(int));
	cudaMemcpy(cuda_params->elem_type_d, elem_type, 
		   nelem * sizeof(int), cudaMemcpyHostToDevice);
}

template<>
void micropp<3>::cuda_finalize()
{
	cudaFree(cuda_params.elem_type_d);
}

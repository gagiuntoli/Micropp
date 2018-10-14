#ifndef __MYWRAPPER_H
#define __MYWRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct micro3 micro3;

	micro3 *init3_(int *ngp, const int size[3], const int *micro_type,
	                   const double *micro_params, const material_t *materials);

#ifdef __cplusplus
}
#endif
#endif

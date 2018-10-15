#ifndef __C_WRAPPER_H
#define __C_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

//	typedef struct micropp3 micropp3;

//	micropp3 *init3_(int *ngp, const int size[3], const int *micro_type,
//	                   const double *micro_params, const material_t *materials);

	typedef struct material_t material_t;
	material_t* micropp_C_material_create();
	void micropp_C_material_set(material_t *in, const double E, double nu,
								double Ka, double Sy, int type);
	void micropp_C_material_print(material_t *in);

#ifdef __cplusplus
}
#endif
#endif

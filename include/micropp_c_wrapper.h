#ifndef C_WRAPPER_H
#define C_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

//	typedef struct micropp3 micropp3;

//	micropp3 *init3_(int *ngp, const int size[3], const int *micro_type,
//	                   const double *micro_params, const material_t *materials);

	void micropp_C_material_create();
	void micropp_C_material_set(int num_mat, const double E, double nu,
								double Ka, double Sy, int type);
	void micropp_C_material_print(int num_mat);

	void micropp_C_create3(int ngp, int size[3], int type, double *params);

#ifdef __cplusplus
}
#endif
#endif

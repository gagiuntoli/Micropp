#ifndef MICROPP_C_WRAPPER_H
#define MICROPP_C_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

    void micropp_C_material_create();
    void micropp_C_material_set(int num_mat, const double E, double nu,
                                double Ka, double Sy, int type);
    void micropp_C_material_print(int num_mat);

    void micropp_C_create3(int ngp, int size[3], int type, double *params);
    void micropp_C_destroy3();
    void micropp_C_set_strain3(int gp, double strain[6]);
    void micropp_C_get_stress3(int gp, double stress[6]);
    void micropp_C_get_ctan3(int gp, double ctan[36]);
    void micropp_C_update_vars();
    void micropp_C_homogenize();
    void micropp_C_print_info();

#ifdef __cplusplus
}
#endif
#endif

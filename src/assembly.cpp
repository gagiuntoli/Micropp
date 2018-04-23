#include "micro.h"

//int assembly_res_ell(double *norm, double *strain_mac)
//{
//  int npe = mesh_struct.npe;
//  int dim = mesh_struct.dim;
//  int nn = mesh_struct.nn;
//  int nelm = mesh_struct.nelm;
//  double *res_e = malloc((dim*npe) * sizeof(double));
//
//  for (int i = 0 ; i < (nn*dim) ; i++) res_ell[i] = 0.0;
//
//  for (int e = 0 ; e < nelm ; e++) {
//
//    for (int i = 0 ; i < (npe*dim) ; i++) res_e[i] = 0.0;
//    mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);
//    
//    for (int gp = 0 ; gp < ngp ; gp++) { // integration res_e
//      get_strain(e, gp, strain_gp);
//      get_stress(e, gp, strain_gp, stress_gp);
//      for (int i = 0 ; i < (npe*dim) ; i++)
//	for (int j = 0 ; j < nvoi ; j++)
//	  res_e[i] += struct_bmat[j][i][gp] * stress_gp[j] * struct_wp[gp];
//    }
//    for (int i = 0 ; i < (npe*dim) ; i++) res_ell[elem_index[i]] += res_e[i]; // assembly
//  }
//
//  if (params.fe2_bc == BC_USTRAIN) {
//    for (int d = 0; d < dim ; d++) {
//      for (int n = 0 ; n < (mesh_struct.ny - 2) ; n++) {
//	res_ell[mesh_struct.nods_x0[n]*dim + d] = 0.0;
//	res_ell[mesh_struct.nods_x1[n]*dim + d] = 0.0;
//      }
//      for (int n = 0 ; n < (mesh_struct.nx - 2) ; n++) {
//	res_ell[mesh_struct.nods_y0[n]*dim + d] = 0.0;
//	res_ell[mesh_struct.nods_y1[n]*dim + d] = 0.0;
//      }
//      res_ell[mesh_struct.nod_x0y0*dim + d] = 0.0;
//      res_ell[mesh_struct.nod_x1y0*dim + d] = 0.0;
//      res_ell[mesh_struct.nod_x1y1*dim + d] = 0.0;
//      res_ell[mesh_struct.nod_x0y1*dim + d] = 0.0;
//    }
//  }
//
//  *norm = 0;
//  for (int i = 0 ; i < (nn*dim) ; i++)
//    *norm += res_ell[i] * res_ell[i];
//  *norm = sqrt(*norm);
//
//  return 0;
//}

void ell_assembly_A (ell_matrix &A, Problem &problem)
{
  int npe = problem.npe;
  int dim = problem.dim;
  int nelm = problem.nelm;
  double Ae[8][8], Ae_vec[64];

  ell_set_zero_mat(&jac_ell);
  for (int e = 0 ; e < nelm ; e++) {

    mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);

    get_elemental (e, Ae, problem);
    for (int i=0; i<npe*dim; i++)
      for (int j=0; j<npe*dim; j++)
	Ae_vec[i*8+j] = Ae[i][j];

    ell_add_vals(A, elem_index, npe*dim, elem_index, npe*dim, Ae_vec); // assembly
  }

//  for (int d = 0; d < dim ; d++) {
//    for (int n = 0; n < mesh_struct.ny - 2 ; n++) {
//      ell_set_zero_row (&jac_ell, mesh_struct.nods_x0[n]*dim + d, 1.0);
//      ell_set_zero_col (&jac_ell, mesh_struct.nods_x0[n]*dim + d, 1.0);
//      ell_set_zero_row (&jac_ell, mesh_struct.nods_x1[n]*dim + d, 1.0);
//      ell_set_zero_col (&jac_ell, mesh_struct.nods_x1[n]*dim + d, 1.0);
//    }
//    for (int n = 0; n < mesh_struct.nx - 2 ; n++) {
//      ell_set_zero_row (&jac_ell, mesh_struct.nods_y0[n]*dim + d, 1.0);
//      ell_set_zero_col (&jac_ell, mesh_struct.nods_y0[n]*dim + d, 1.0);
//      ell_set_zero_row (&jac_ell, mesh_struct.nods_y1[n]*dim + d, 1.0);
//      ell_set_zero_col (&jac_ell, mesh_struct.nods_y1[n]*dim + d, 1.0);
//    }
//    ell_set_zero_row (&jac_ell, mesh_struct.nod_x0y0*dim + d, 1.0);
//    ell_set_zero_col (&jac_ell, mesh_struct.nod_x0y0*dim + d, 1.0);
//    ell_set_zero_row (&jac_ell, mesh_struct.nod_x1y0*dim + d, 1.0);
//    ell_set_zero_col (&jac_ell, mesh_struct.nod_x1y0*dim + d, 1.0);
//    ell_set_zero_row (&jac_ell, mesh_struct.nod_x1y1*dim + d, 1.0);
//    ell_set_zero_col (&jac_ell, mesh_struct.nod_x1y1*dim + d, 1.0);
//    ell_set_zero_row (&jac_ell, mesh_struct.nod_x0y1*dim + d, 1.0);
//    ell_set_zero_col (&jac_ell, mesh_struct.nod_x0y1*dim + d, 1.0);
//  }

}

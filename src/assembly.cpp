#include <cmath>
#include "micro.h"

double assembly_b (double *b, Problem &problem)
{
  int npe = problem.npe;
  int dim = problem.dim;
  int nx = problem.nx;
  int ny = problem.ny;
  int nn = problem.nn;
  int nelem = problem.nelem;
  int index[8];
  double be[8];

  for (int i=0 ; i<nn*dim; i++) {
    b[i] = 0.0;
  }

  for (int e=0; e<nelem; e++) {

    int xfactor = e%(nx-1);
    int yfactor = e/(ny-1);

    int n0 = yfactor     * nx + xfactor     ;
    int n1 = yfactor     * nx + xfactor + 1 ;
    int n2 = (yfactor+1) * nx + xfactor + 1 ;
    int n3 = (yfactor+1) * nx + xfactor     ;
    index[0] = n0*dim; index[1] = n0*dim + 1;
    index[2] = n1*dim; index[3] = n1*dim + 1;
    index[4] = n2*dim; index[5] = n2*dim + 1;
    index[6] = n3*dim; index[7] = n3*dim + 1;
    //cout << "e : n0="<<n0<<" n1="<<n1<<" n2="<<n2<<" n3="<<n3<<endl;
    
    get_elemental_b(e, be, problem);
         
//    for (int gp = 0 ; gp < ngp ; gp++) { // integration res_e
//      get_strain(e, gp, strain_gp);
//      get_stress(e, gp, strain_gp, stress_gp);
//      for (int i = 0 ; i < (npe*dim) ; i++)
//	for (int j = 0 ; j < nvoi ; j++)
//	  res_e[i] += struct_bmat[j][i][gp] * stress_gp[j] * struct_wp[gp];
//    }
    for (int i=0 ; i<npe*dim; i++) {
      b[index[i]] += be[i]; // assembly
    }
  }

//  for (int d = 0; d < dim ; d++) {
//    for (int n = 0 ; n < (mesh_struct.ny - 2) ; n++) {
//      b[problem.nods_x0[n]*dim + d] = 0.0;
//      b[problem.nods_x1[n]*dim + d] = 0.0;
//    }
//    for (int n = 0 ; n < (problem.nx - 2) ; n++) {
//      b[problem.nods_y0[n]*dim + d] = 0.0;
//      b[problem.nods_y1[n]*dim + d] = 0.0;
//    }
//    b[problem.nod_x0y0*dim + d] = 0.0;
//    b[problem.nod_x1y0*dim + d] = 0.0;
//    b[problem.nod_x1y1*dim + d] = 0.0;
//    b[problem.nod_x0y1*dim + d] = 0.0;
//  }

  double norm = 0.0;
  for (int i=0; i<nn*dim; i++) {
    norm += b[i]*b[i];
  }
  norm = sqrt(norm);

  return norm;
}

void assembly_A (ell_matrix &A, Problem &problem)
{
  int npe = problem.npe;
  int dim = problem.dim;
  int nx = problem.nx;
  int ny = problem.ny;
  int nelem = problem.nelem;
  int index[8];
  double Ae[8][8], Ae_vec[64];

  ell_set_zero_mat(&A);
  for (int e = 0 ; e < nelem ; e++) {

    int xfactor = e%(nx-1);
    int yfactor = e/(ny-1);

    int n0 = yfactor     * nx + xfactor     ;
    int n1 = yfactor     * nx + xfactor + 1 ;
    int n2 = (yfactor+1) * nx + xfactor + 1 ;
    int n3 = (yfactor+1) * nx + xfactor     ;
    index[0] = n0*dim; index[1] = n0*dim + 1;
    index[2] = n1*dim; index[3] = n1*dim + 1;
    index[4] = n2*dim; index[5] = n2*dim + 1;
    index[6] = n3*dim; index[7] = n3*dim + 1;
    //cout << "e : n0="<<n0<<" n1="<<n1<<" n2="<<n2<<" n3="<<n3<<endl;

    get_elemental (e, Ae, problem);
    for (int i=0; i<npe*dim; i++)
      for (int j=0; j<npe*dim; j++)
	Ae_vec[i*8+j] = Ae[i][j];

    ell_add_vals(&A, index, npe*dim, index, npe*dim, Ae_vec); // assembly
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

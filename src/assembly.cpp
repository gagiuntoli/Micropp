#include <cmath>
#include <iostream>
#include "micro.h"

using namespace std;

void Problem::setDisp (double *eps)
{
  /* dim 2 
     | eps(0)    eps(2)/2 |
     | eps(2)/2  eps(1)   |
  */
  double dux, duy, xcoor, ycoor;
  
  // y = 0
  for (int i=0; i<nx; i++) {
    xcoor = i*dx;
    ycoor = 0.0;
    dux =     eps[0]*xcoor + 0.5*eps[2]*ycoor;
    duy = 0.5*eps[2]*xcoor +     eps[1]*ycoor;
    u[i*dim  ] = dux;
    u[i*dim+1] = duy;
  }
  // y = ly
  for (int i=0; i<nx; i++) {
    xcoor = i*dx;
    ycoor = ly;
    dux =     eps[0]*xcoor + 0.5*eps[2]*ycoor;
    duy = 0.5*eps[2]*xcoor +     eps[1]*ycoor;
    u[(i+(ny-1)*nx)*dim  ] = dux;
    u[(i+(ny-1)*nx)*dim+1] = duy;
  }
  // x = 0
  for (int i=0; i<ny-2; i++) {
    xcoor = 0.0;
    ycoor = (i+1)*dy;
    dux =     eps[0]*xcoor + 0.5*eps[2]*ycoor;
    duy = 0.5*eps[2]*xcoor +     eps[1]*ycoor;
    u[(i+1)*nx*dim  ] = dux;
    u[(i+1)*nx*dim+1] = duy;
  }
  // x = lx
  for (int i=0; i<ny-2; i++) {
    xcoor = lx;
    ycoor = (i+1)*dy;
    dux =       eps[0]*xcoor + 0.5 * eps[2]*ycoor;
    duy = 0.5 * eps[2]*xcoor +       eps[1]*ycoor;
    u[((i+2)*nx-1)*dim  ] = dux;
    u[((i+2)*nx-1)*dim+1] = duy;
  }

  if (flag_print_u == true) {
    for (int i=0; i<nn; i++) {
      cout << u[i*dim] << " " << u[i*dim+1] << endl;
    }
  }

}

double Problem::Assembly_b (void)
{
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
    
    getElemental_b (e, be);
         
    for (int i=0; i<npe*dim; i++) {
      b[index[i]] += be[i]; // assembly
    }
  } // loop over elements

  // boundary conditions
  // y = 0
  for (int i=0; i<nx; i++) {
    for (int d=0; d<dim; d++) {
      b[i*dim + d] = 0.0;
    }
  }
  // y = ly
  for (int i=0; i<nx; i++) {
    for (int d=0; d<dim; d++) {
      b[(i+(ny-1)*nx)*dim + d] = 0.0;
    }
  }
  // x = 0
  for (int i=0; i<ny-2; i++) {
    for (int d=0; d<dim; d++) {
      b[(i+1)*nx*dim + d] = 0.0;
    }
  }
  // x = lx
  for (int i=0; i<ny-2; i++) {
    for (int d=0; d<dim; d++) {
      b[((i+2)*nx-1)*dim + d] = 0.0;
    }
  }
  // end of boundary conditions

  if (flag_print_b == true) {
    for (int i=0; i<nn*dim; i++) {
      cout << b[i] << endl;
    }
  }

  double norm = 0.0;
  for (int i=0; i<nn*dim; i++) {
    norm += b[i]*b[i];
  }
  norm = sqrt(norm);

  return norm;
}

void Problem::Assembly_A (void)
{
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

    get_elemental_A (e, Ae);
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
  if (flag_print_A == true) {
    ell_print (&A);
  }

}

void Problem::get_elemental_A (int e, double (&Ae)[8][8])
{
/*   // for debugging
     for (int i=0; i<8; i++)
     for (int j=0; j<8; j++)
     Ae[i][j] = 1;
 */
  
  double nu = 0.3, E;
  double ctan[3][3];
  E  = 1e7;
  if (distance(e) < -1.75) {
    E  = 1.0e7;
  }
  else {
     E  = 1.0e6;
  }
  ctan[0][0]=(1-nu); ctan[0][1]=nu    ; ctan[0][2]=0;
  ctan[1][0]=nu    ; ctan[1][1]=(1-nu); ctan[1][2]=0;
  ctan[2][0]=0     ; ctan[2][1]=0     ; ctan[2][2]=(1-2*nu)/2;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      ctan[i][j] *= E/((1+nu)*(1-2*nu));
    }
  }

  double xg[4][2] = {
    {-0.577350269189626, -0.577350269189626},
    {+0.577350269189626, -0.577350269189626},
    {+0.577350269189626, +0.577350269189626},
    {-0.577350269189626, +0.577350269189626}};

  double dsh[4][2], b_mat[3][8], cxb[3][8];

  for (int i=0; i<8; i++) 
    for (int j=0; j<8; j++) 
      Ae[i][j] = 0.0;

  for (int gp=0; gp<4; gp++) {

    dsh[0][0] = -(1-xg[gp][1])/4*2/dx ; dsh[0][1] = -(1-xg[gp][0])/4*2/dy;
    dsh[1][0] = +(1-xg[gp][1])/4*2/dx ; dsh[1][1] = -(1+xg[gp][0])/4*2/dy;
    dsh[2][0] = +(1+xg[gp][1])/4*2/dx ; dsh[2][1] = +(1+xg[gp][0])/4*2/dy;
    dsh[3][0] = -(1+xg[gp][1])/4*2/dx ; dsh[3][1] = +(1-xg[gp][0])/4*2/dy;

    for (int i=0; i<4; i++) {
      b_mat[1][i*dim] = dsh[i][0]; b_mat[1][i*dim + 1] = 0        ;
      b_mat[2][i*dim] = 0        ; b_mat[2][i*dim + 1] = dsh[i][1];
      b_mat[3][i*dim] = dsh[i][1]; b_mat[3][i*dim + 1] = dsh[i][0];
    }

    for (int i=0; i<3; i++) {
      for (int j=0; j<8; j++) {
	for (int k=0; k<3; k++) {
	  cxb[i][j] += ctan[i][k] * b_mat[k][j];
	}
      }
    }

    for (int i=0; i<8; i++) {
      for (int j=0; j<8; j++) {
	for (int m=0; m<3; m++) {
	    Ae[i][j] += b_mat[m][i] * cxb[m][j] * wg[gp];
	}
      }
    }
  } // gp loop
}

double Problem::distance (int e)
{
  int xfactor = e%(nx-1);
  int yfactor = e/(ny-1);
  double xdis = pow(xfactor*dx + dx/2 - lx/2,2);
  double ydis = pow(yfactor*dy + dy/2 - ly/2,2);
  return sqrt(xdis + ydis);
}

void Problem::getElemental_b (int e, double (&be)[8])
{
  /*
     // for debugging
     for (int i=0; i<8; i++)
        be[i] = 1;
   */

  double dsh[4][2], bmat[3][8], cxb[3][8], stress_gp[6];
  double xg[4][2] = {
    {-0.577350269189626, -0.577350269189626},
    {+0.577350269189626, -0.577350269189626},
    {+0.577350269189626, +0.577350269189626},
    {-0.577350269189626, +0.577350269189626}};

  for (int i=0; i<8; i++) 
    be[i] = 0.0;

  for (int gp=0; gp<4; gp++) {

    dsh[0][0] = -(1-xg[gp][1])/4*2/dx ; dsh[0][1] = -(1-xg[gp][0])/4*2/dy;
    dsh[1][0] = +(1-xg[gp][1])/4*2/dx ; dsh[1][1] = -(1+xg[gp][0])/4*2/dy;
    dsh[2][0] = +(1+xg[gp][1])/4*2/dx ; dsh[2][1] = +(1+xg[gp][0])/4*2/dy;
    dsh[3][0] = -(1+xg[gp][1])/4*2/dx ; dsh[3][1] = +(1-xg[gp][0])/4*2/dy;

    for (int i=0; i<4; i++) {
      bmat[0][i*dim] = dsh[i][0]; bmat[0][i*dim+1] = 0        ;
      bmat[1][i*dim] = 0        ; bmat[1][i*dim+1] = dsh[i][1];
      bmat[2][i*dim] = dsh[i][1]; bmat[2][i*dim+1] = dsh[i][0];
    }

    getStress (e, gp, stress_gp);

    for (int i=0; i<npe*dim; i++) {
      for (int j=0; j<nvoi; j++) {
	be[i] += bmat[j][i] * stress_gp[j] * wg[gp];
      }
    }

  } // gp loop
}

void Problem::getStress (int e, int gp, double *stress_gp)
{
  if (dim == 2) {
    double strain_gp[3];
    getStrain (e, gp, strain_gp);

    double nu = 0.3, E;
    double ctan[3][3];

    if (elem_type[e] == 0) {
      E  = 1.0e6;
    } else {
      E  = 1.0e7;
    }
    ctan[0][0]=(1-nu); ctan[0][1]=nu    ; ctan[0][2]=0;
    ctan[1][0]=nu    ; ctan[1][1]=(1-nu); ctan[1][2]=0;
    ctan[2][0]=0     ; ctan[2][1]=0     ; ctan[2][2]=(1-2*nu)/2;
    for (int i=0; i<nvoi; i++) {
      for (int j=0; j<nvoi; j++) {
	ctan[i][j] *= E/((1+nu)*(1-2*nu));
      }
    }

    for (int i=0; i<nvoi; i++) { 
      for (int j=0; j<nvoi; j++) { 
	stress_gp[i] = ctan[i][j] * strain_gp[j];
      }
    }

  }
}

void Problem::getStrain (int e, int gp, double *strain_gp)
{
  double elem_disp[3*8];
  getElemDisp (e, elem_disp);

  if (dim == 2) {

    double xg[4][2] = {
      {-0.577350269189626, -0.577350269189626},
      {+0.577350269189626, -0.577350269189626},
      {+0.577350269189626, +0.577350269189626},
      {-0.577350269189626, +0.577350269189626}};

    double dsh[4][2], bmat[3][8], cxb[3][8];

    dsh[0][0] = -(1-xg[gp][1])/4*2/dx ; dsh[0][1] = -(1-xg[gp][0])/4*2/dy;
    dsh[1][0] = +(1-xg[gp][1])/4*2/dx ; dsh[1][1] = -(1+xg[gp][0])/4*2/dy;
    dsh[2][0] = +(1+xg[gp][1])/4*2/dx ; dsh[2][1] = +(1+xg[gp][0])/4*2/dy;
    dsh[3][0] = -(1+xg[gp][1])/4*2/dx ; dsh[3][1] = +(1-xg[gp][0])/4*2/dy;

    for (int i=0; i<4; i++) {
      bmat[0][i*dim] = dsh[i][0]; bmat[0][i*dim + 1] = 0        ;
      bmat[1][i*dim] = 0        ; bmat[1][i*dim + 1] = dsh[i][1];
      bmat[2][i*dim] = dsh[i][1]; bmat[2][i*dim + 1] = dsh[i][0];
    }
    for (int v=0; v<nvoi; v++) {
      strain_gp[v] = 0.0;
      for (int i=0; i<npe*dim; i++) {
	strain_gp[v] += bmat[v][i] * elem_disp[i];
      }
    }

  }

}

void Problem::getElemDisp (int e, double *elem_disp)
{
  if (dim == 2) {
    int xfactor = e%(nx-1);
    int yfactor = e/(ny-1);

    int n0 = yfactor     * nx + xfactor     ;
    int n1 = yfactor     * nx + xfactor + 1 ;
    int n2 = (yfactor+1) * nx + xfactor + 1 ;
    int n3 = (yfactor+1) * nx + xfactor     ;

    for (int d=0; d<dim; d++) {
      elem_disp[0*dim + d] = u[n0*dim + d];
      elem_disp[1*dim + d] = u[n1*dim + d];
      elem_disp[2*dim + d] = u[n2*dim + d];
      elem_disp[3*dim + d] = u[n3*dim + d];
    }
  }
}

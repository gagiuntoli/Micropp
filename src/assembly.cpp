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
  
  // y = 0
  for (int i=0; i<nx; i++) {
    double xcoor = i*dx;
    double ycoor = 0.0;
    double dux =     eps[0]*xcoor + 0.5*eps[2]*ycoor;
    double duy = 0.5*eps[2]*xcoor +     eps[1]*ycoor;
    u[i*dim  ] = dux;
    u[i*dim+1] = duy;
  }
  // y = ly
  for (int i=0; i<nx; i++) {
    double xcoor = i*dx;
    double ycoor = ly;
    double dux =     eps[0]*xcoor + 0.5*eps[2]*ycoor;
    double duy = 0.5*eps[2]*xcoor +     eps[1]*ycoor;
    u[(i+(ny-1)*nx)*dim  ] = dux;
    u[(i+(ny-1)*nx)*dim+1] = duy;
  }
  // x = 0
  for (int i=0; i<ny-2; i++) {
    double xcoor = 0.0;
    double ycoor = (i+1)*dy;
    double dux =     eps[0]*xcoor + 0.5*eps[2]*ycoor;
    double duy = 0.5*eps[2]*xcoor +     eps[1]*ycoor;
    u[(i+1)*nx*dim  ] = dux;
    u[(i+1)*nx*dim+1] = duy;
  }
  // x = lx
  for (int i=0; i<ny-2; i++) {
    double xcoor = lx;
    double ycoor = (i+1)*dy;
    double dux =       eps[0]*xcoor + 0.5 * eps[2]*ycoor;
    double duy = 0.5 * eps[2]*xcoor +       eps[1]*ycoor;
    u[((i+2)*nx-1)*dim  ] = dux;
    u[((i+2)*nx-1)*dim+1] = duy;
  }

  if (flag_print_u == true)
    for (int i=0; i<nn; i++)
      cout << u[i*dim] << " " << u[i*dim+1] << endl;
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
  for (int i=0; i<nx; i++)
    for (int d=0; d<dim; d++)
      b[i*dim + d] = 0.0;

  // y = ly
  for (int i=0; i<nx; i++)
    for (int d=0; d<dim; d++)
      b[(i+(ny-1)*nx)*dim + d] = 0.0;

  // x = 0
  for (int i=0; i<ny-2; i++)
    for (int d=0; d<dim; d++)
      b[(i+1)*nx*dim + d] = 0.0;

  // x = lx
  for (int i=0; i<ny-2; i++)
    for (int d=0; d<dim; d++)
      b[((i+2)*nx-1)*dim + d] = 0.0;

  // end of boundary conditions

  if (flag_print_b == true)
    for (int i=0; i<nn*dim; i++)
      cout << b[i*dim] << " " << b[i*dim+1] << endl;

  for (int i=0 ; i<nn*dim; i++) {
    b[i] = -b[i];
  }

  double norm = 0.0;
  for (int i=0; i<nn*dim; i++)
    norm += b[i]*b[i];
  norm = sqrt(norm);

  return norm;
}

void Problem::Assembly_A (void)
{
  int index[8];
  double Ae[64];

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

    getElemental_A (e, Ae);

//    ell_add_vals(&A, index, npe*dim, index, npe*dim, Ae); // assembly
    ell_add_2D (A, e, Ae, dim, nx, ny);
  }

  // y = 0
  for (int i=0; i<nx; i++) {
    for (int d=0; d<dim; d++) {
      int index = i*dim + d;
      ell_set_zero_row (&A, index, 1.0);
      ell_set_zero_col (&A, index, 1.0);
    }
  }
  // y = ly
  for (int i=0; i<nx; i++) {
    for (int d=0; d<dim; d++) {
      int index = (i+(ny-1)*nx)*dim + d;
      ell_set_zero_row (&A, index, 1.0);
      ell_set_zero_col (&A, index, 1.0);
    }
  }
  // x = 0
  for (int i=0; i<ny-2; i++) {
    for (int d=0; d<dim; d++) {
      int index = (i+1)*nx*dim + d;
      ell_set_zero_row (&A, index, 1.0);
      ell_set_zero_col (&A, index, 1.0);
    }
  }
  // x = lx
  for (int i=0; i<ny-2; i++) {
    for (int d=0; d<dim; d++) {
      int index = ((i+2)*nx-1)*dim + d;
      ell_set_zero_row (&A, index, 1.0);
      ell_set_zero_col (&A, index, 1.0);
    }
  }

  if (flag_print_A == true) {
    ell_print (&A);
  }

}

void Problem::getElemental_A (int e, double (&Ae)[64])
{
/*   // for debugging
     for (int i=0; i<8; i++)
     for (int j=0; j<8; j++)
     Ae[i][j] = 1;
 */
  
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

  for (int i=0; i<nvoi; i++)
    for (int j=0; j<nvoi; j++)
      ctan[i][j] *= E/((1+nu)*(1-2*nu));

  double xg[4][2] = {
    {-0.577350269189626, -0.577350269189626},
    {+0.577350269189626, -0.577350269189626},
    {+0.577350269189626, +0.577350269189626},
    {-0.577350269189626, +0.577350269189626}};

  double dsh[4][2], b_mat[3][8], cxb[3][8];

  for (int i=0; i<npe*dim*npe*dim; i++) 
    Ae[i] = 0.0;

  for (int gp=0; gp<4; gp++) {

    dsh[0][0] = -(1-xg[gp][1])/4*2/dx ; dsh[0][1] = -(1-xg[gp][0])/4*2/dy;
    dsh[1][0] = +(1-xg[gp][1])/4*2/dx ; dsh[1][1] = -(1+xg[gp][0])/4*2/dy;
    dsh[2][0] = +(1+xg[gp][1])/4*2/dx ; dsh[2][1] = +(1+xg[gp][0])/4*2/dy;
    dsh[3][0] = -(1+xg[gp][1])/4*2/dx ; dsh[3][1] = +(1-xg[gp][0])/4*2/dy;

    for (int i=0; i<4; i++) {
      b_mat[0][i*dim] = dsh[i][0]; b_mat[0][i*dim + 1] = 0        ;
      b_mat[1][i*dim] = 0        ; b_mat[1][i*dim + 1] = dsh[i][1];
      b_mat[2][i*dim] = dsh[i][1]; b_mat[2][i*dim + 1] = dsh[i][0];
    }

    for (int i=0; i<nvoi; i++) {
      for (int j=0; j<npe*dim; j++) {
	cxb[i][j] = 0.0;
	for (int k=0; k<nvoi; k++)
	  cxb[i][j] += ctan[i][k] * b_mat[k][j];
      }
    }

    double wg = 0.25*dx*dy;
    for (int i=0; i<npe*dim; i++)
      for (int j=0; j<npe*dim; j++)
	for (int m=0; m<nvoi; m++)
	    Ae[i*npe*dim + j] += b_mat[m][i] * cxb[m][j] * wg;

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

    double wg = 0.25*dx*dy;
    for (int i=0; i<npe*dim; i++) {
      for (int j=0; j<nvoi; j++) {
	be[i] += bmat[j][i] * stress_gp[j] * wg;
      }
    }

  } // gp loop
}

void Problem::calcAverageStress (void)
{
  for (int v=0; v<nvoi; v++)
    stress_ave[v] = 0.0;

  for (int e=0; e<nelem; e++) {

    double stress_aux[6];
    for (int i=0; i<nvoi; i++)
      stress_aux[i] = 0.0;

    for (int gp=0; gp<4; gp++) {

      double stress_gp[6];
      double wg = 0.25*dx*dy;

      getStress (e, gp, stress_gp);
      for (int v=0; v<nvoi; v++)
	stress_aux[v] += stress_gp[v] * wg;

    }
    for (int v=0; v<nvoi; v++)
      stress_ave[v] += stress_aux[v];
  }

  for (int v=0; v<nvoi; v++)
    stress_ave[v] /= (lx*ly);
}

void Problem::calcAverageStrain (void)
{
  for (int v=0; v<nvoi; v++)
    strain_ave[v] = 0.0;

  for (int e=0; e<nelem; e++) {

    double strain_aux[6];
    for (int i=0; i<nvoi; i++)
      strain_aux[i] = 0.0;

    for (int gp=0; gp<4; gp++) {

      double strain_gp[6];
      double wg = 0.25*dx*dy;

      getStrain (e, gp, strain_gp);
      for (int v=0; v<nvoi; v++)
	strain_aux[v] += strain_gp[v] * wg;

    }

    for (int v=0; v<nvoi; v++)
      strain_ave[v] += strain_aux[v];
  }

  for (int v=0; v<nvoi; v++)
    strain_ave[v] /= (lx*ly);
}

void Problem::calcDistributions (void)
{
  for (int e=0; e<nelem; e++) {

    double strain_aux[6];
    double stress_aux[6];
    for (int i=0; i<nvoi; i++) {
      strain_aux[i] = 0.0;
      stress_aux[i] = 0.0;
    }

    for (int gp=0; gp<4; gp++) {

      double stress_gp[6], strain_gp[6];
      double wg = 0.25*dx*dy;

      getStress (e, gp, stress_gp);
      getStrain (e, gp, strain_gp);
      for (int v=0; v<nvoi; v++) {
	strain_aux[v] += strain_gp[v] * wg;
	stress_aux[v] += stress_gp[v] * wg;
      }

    }
    double vol = dx*dy;
    for (int v=0; v<nvoi; v++) {
      strain[e*nvoi + v] = strain_aux[v] / vol;
      stress[e*nvoi + v] = stress_aux[v] / vol;
    }
  }
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
    for (int i=0; i<nvoi; i++)
      for (int j=0; j<nvoi; j++)
	ctan[i][j] *= E/((1+nu)*(1-2*nu));

    for (int i=0; i<nvoi; i++) {
      stress_gp[i] = 0.0;
      for (int j=0; j<nvoi; j++) {
	stress_gp[i] += ctan[i][j] * strain_gp[j];
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

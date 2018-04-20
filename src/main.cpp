#include <iostream>
#include <ctime>
#include <cmath>
#include "micro.h"
#include "csr.h"

using namespace std;

int csr_set_A (csr_matrix &A, Problem &problem)
{
  // set offsets and cols
  int dim = problem.dim;
  int npe = problem.npe;
  int nelem = problem.nelem; 
  int nx = problem.nx; 
  int ny = problem.ny; 
  int nn = problem.nn; 

  A.row_offsets[0] = 0;
  for (int i=1; i<(A.num_rows+1); i++)
    A.row_offsets[i] =  18*i;
  
  for (int i=0; i<nn; i++){
    // the coorners
    if (i == 0) {       
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
	}
      }
    }
    else if (i == nx-1) {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
	}
      }
    }
    else if (i == nx*ny-1) {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
	}
      }
    }
    else if (i == (ny-1)*nx) {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
	}
      }
    }
    // y=0
    else if (i < nx) {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
	}
      }
    }
    // y=ly
    else if (i > (ny-1)*nx) {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
	}
      }
    }
    // x=0
    else if ((i%nx) == 0) {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
	}
      }
    }
    // x=ly
    else if ((i+1)%nx == 0) {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0; 
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = 0;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
	}
      }
    }
    // not y=0 neither y=ly, x=0, x=lx
    else {
      for (int d1=0; d1<2; d1++){
	for (int d=0; d<2; d++){
	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
	}
      }
    }
  }

//  cout << "A.num_rows = " << A.num_rows << endl;
//  for (int i=0; i<(A.num_rows+1); i++)
//    cout << A.row_offsets[i] << " ";
//  cout << endl << endl;
//  for (int i=0; i<A.num_rows; i++) {
//    for (int j=0; j<18; j++)
//      cout << A.cols[i*18 + j] << " ";
//    cout << endl;
//  }
//  cout << endl ;

}

double distance (int e, Problem &problem)
{
  int nx = problem.nx;
  int ny = problem.ny;
  double dx = problem.dx;
  double dy = problem.dy;
  double lx = problem.lx;
  double ly = problem.ly;
  int xfactor = e%(nx-1);
  int yfactor = e/(ny-1);
  double xdis = pow(xfactor*dx + dx/2 - lx/2,2);
  double ydis = pow(yfactor*dy + dy/2 - ly/2,2);
  return sqrt(xdis + ydis);
}

void get_elemental (int e, double (&Ae)[8][8], Problem &problem)
{
  /* // for debugging purposes
     for (int i=0; i<8; i++)
     for (int j=0; j<8; j++)
     Ae[i][j] = 1;
   */

  double nu = 0.3, E;
  double ctan[3][3];
  E  = 1e7;
  if (distance(e, problem) < -1.75) {
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

  for (int i=0; i<8; i++) {
    for (int j=0; j<8; j++) {
      Ae[i][j] = 0.0;
      for (int gp=0; gp<4; gp++) {
	for (int m=0; m<8; m++) {
	  for (int k=0; k<3; k++) {
	    Ae[i][j] += problem.b_mat[m][i][gp] * ctan[m][k] * problem.b_mat[k][j][gp] * problem.wg[gp];
	  }
	}
      }
    }
  }
  
}

void csr_assembly_A (csr_matrix &A, Problem &problem)
{
  int dim = problem.dim;
  int npe = problem.npe;
  int nelem = problem.nelem; 
  int nx = problem.nx; 
  int ny = problem.ny; 

  int n0, n1, n2, n3;
  int cols_row_0[4]={4,5,8,7};
  int cols_row_1[4]={3,4,7,6};
  int cols_row_2[4]={0,1,4,3};
  int cols_row_3[4]={1,2,5,4};
  int yfactor;
  int xfactor;
  int index;

  double Ae[8][8];
  clock_t start, end, start_1, end_1;
  double elapsed_secs, elaps_1=0.0, elaps_2=0.0;

  for (int i=0; i<A.num_rows*18; i++)
    A.coefs[i] = 0.0;

  start = clock();
  for (int e=0; e<nelem; e++){

    xfactor = e%(nx-1);
    yfactor = e/(ny-1);

    n0 = yfactor     * nx + xfactor     ;
    n1 = yfactor     * nx + xfactor + 1 ;
    n2 = (yfactor+1) * nx + xfactor + 1 ;
    n3 = (yfactor+1) * nx + xfactor     ;
    //cout << "e : n0="<<n0<<" n1="<<n1<<" n2="<<n2<<" n3="<<n3<<endl;

    start_1 = clock();
    get_elemental (e, Ae, problem);
    end_1 = clock();
    elaps_1 += double(end_1 - start_1) / CLOCKS_PER_SEC;

    start_1 = clock();
    for (int n=0; n<4; n++){
      for (int d=0; d<2; d++){
	A.coefs[n0*dim*18 + cols_row_0[n]*dim + d     ] += Ae[0][n*2+d];
	A.coefs[n0*dim*18 + cols_row_0[n]*dim + d + 18] += Ae[1][n*2+d];
	A.coefs[n1*dim*18 + cols_row_1[n]*dim + d     ] += Ae[2][n*2+d];
	A.coefs[n1*dim*18 + cols_row_1[n]*dim + d + 18] += Ae[3][n*2+d];
	A.coefs[n2*dim*18 + cols_row_2[n]*dim + d     ] += Ae[4][n*2+d];
	A.coefs[n2*dim*18 + cols_row_2[n]*dim + d + 18] += Ae[5][n*2+d];
	A.coefs[n3*dim*18 + cols_row_3[n]*dim + d     ] += Ae[6][n*2+d];
	A.coefs[n3*dim*18 + cols_row_3[n]*dim + d + 18] += Ae[7][n*2+d];
      }
    }
    end_1 = clock();
    elaps_2 += double(end_1 - start_1) / CLOCKS_PER_SEC;

  }
  end = clock();
  elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
  cout << "time assembly loop in e   = "<<elapsed_secs<<endl;
  cout << "time assembly calc     Ae = "<<elaps_1<<endl;
  cout << "time assembly assembly Ae = "<<elaps_2<<endl;

  /* boundary conditions */
  start = clock();
  for (int n=0; n<nx; n++){ // y=0
    for (int k=0; k<18; k++){
      A.coefs[n*dim*18 + k     ] = 0.0;
      A.coefs[n*dim*18 + k + 18] = 0.0;
    }
    A.coefs[n*dim*18 + 8     ] = 1.0;
    A.coefs[n*dim*18 + 9 + 18] = 1.0;
  }
  for (int n=(ny-1)*nx; n<nx*ny; n++){ // y=ly
    for (int k=0; k<18; k++){
      A.coefs[n*dim*18 + k     ] = 0.0;
      A.coefs[n*dim*18 + k + 18] = 0.0;
    }
    A.coefs[n*dim*18 + 8     ] = 1.0;
    A.coefs[n*dim*18 + 9 + 18] = 1.0;
  }
  for (int n=1; n<ny-1; n++){ // x=0
    for (int k=0; k<18; k++){
      A.coefs[(n*nx)*dim*18 + k     ] = 0.0;
      A.coefs[(n*nx)*dim*18 + k + 18] = 0.0;
    }
    A.coefs[(n*nx)*dim*18 + 8     ] = 1.0;
    A.coefs[(n*nx)*dim*18 + 9 + 18] = 1.0;
  }
  for (int n=2; n<ny; n++){ // x=lx
    for (int k=0; k<18; k++){
      A.coefs[(n*nx-1)*dim*18 + k     ] = 0.0;
      A.coefs[(n*nx-1)*dim*18 + k + 18] = 0.0;
    }
    A.coefs[(n*nx-1)*dim*18 + 8     ] = 1.0;
    A.coefs[(n*nx-1)*dim*18 + 9 + 18] = 1.0;
  }
  end = clock();
  elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
  cout << "time assembly bc          = "<<elapsed_secs<<endl;

//  for (int i=0; i<A.num_rows; i++) {
//    for (int j=0; j<18; j++)
//      cout << A.coefs[i*18 + j] << " ";
//    cout << endl;
//  }
//  cout << endl ;

}

void csr_assembly_res (csr_vector &res, Problem &problem)
{
  int dim = problem.dim;
  int nn = problem.nn; 
  for (int i=0; i<nn*dim; i++)
    res.coefs[i] = 1.1;
}

void csr_assembly_dx (csr_vector &dx, Problem &problem)
{
  int dim = problem.dim;
  int nn = problem.nn; 
  for (int i=0; i<nn*dim; i++)
    dx.coefs[i] = 0.0;
}

int main (int argc, char *argv[])
{
  try {

    csr_matrix A;
    csr_vector x, dx, res;
    Problem problem (argc, argv);
    clock_t start, end;
    double elapsed_secs;

    double eps [] = { 0.005, 0.0, 0.0 };
    int nn = problem.nn;
    int dim = problem.dim;

    csr_alloc_A (A, nn*dim);
    csr_alloc_v (res, nn*dim);
    csr_alloc_v (dx, nn*dim);
    csr_alloc_v (x, nn*dim);
    csr_set_A (A, problem);

    cout << "Init micropp -> preparing ..." << endl;
    cout << "dim   = " << problem.dim << endl;
    cout << "nx    = " << problem.nx << endl;
    cout << "ny    = " << problem.ny << endl;
    cout << "nn    = " << problem.nn << endl;
    cout << "nelem = " << problem.nelem << endl;
    cout << "dx    = " << problem.dx << endl;
    cout << "dy    = " << problem.dy << endl;

    // assembly
    start = clock();

    csr_assembly_A   (A, problem);
    csr_assembly_res (res, problem);
    csr_assembly_res (dx, problem);

    end = clock();

    elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    cout << "time assembly             = "<<elapsed_secs<<endl;

    // solve system  A * dx = -res
    start = clock();
    csr_cg (A, dx, res);
    end = clock();
    elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    cout << "time solve                = "<<elapsed_secs<<endl;

    // free memory
    csr_free_A (A);
    csr_free_v (res);
    csr_free_v (dx);
    csr_free_v (x);

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}

#include <../inc/micro.hpp>
#include <../inc/csr.hpp>
#include <iostream>

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
  
  for (int i=0; i<problem.nn; i++){
    if (i>=nx && i<(ny-1)*nx) { // not y=0 neither y=ly
      if ((i%nx)!=0 && (i%(nx-1))!=0) { // not x=0 neither x=lx
	for (int d=0; d<2; d++){
	  A.cols[i*dim*9 + 0*dim + d] = (i - nx - 1)*dim + d;
	  A.cols[i*dim*9 + 1*dim + d] = (i - nx    )*dim + d;
	  A.cols[i*dim*9 + 2*dim + d] = (i - nx + 1)*dim + d;
	  A.cols[i*dim*9 + 3*dim + d] = (i - 1     )*dim + d;
	  A.cols[i*dim*9 + 4*dim + d] = (i         )*dim + d;
	  A.cols[i*dim*9 + 5*dim + d] = (i + 1     )*dim + d;
	  A.cols[i*dim*9 + 6*dim + d] = (i + nx - 1)*dim + d;
	  A.cols[i*dim*9 + 7*dim + d] = (i + nx    )*dim + d;
	  A.cols[i*dim*9 + 8*dim + d] = (i + nx + 1)*dim + d;
	}
      }
    }
    else if (i==0) {       // the coorners
      for (int d=0; d<2; d++){
	A.cols[i*dim*9 + 0*dim + d] = 0;
	A.cols[i*dim*9 + 1*dim + d] = 0;
	A.cols[i*dim*9 + 2*dim + d] = 0;
	A.cols[i*dim*9 + 3*dim + d] = 0;
	A.cols[i*dim*9 + 4*dim + d] = (i         )*dim + d;
	A.cols[i*dim*9 + 5*dim + d] = (i + 1     )*dim + d;
	A.cols[i*dim*9 + 6*dim + d] = 0;
	A.cols[i*dim*9 + 7*dim + d] = (i + nx    )*dim + d;
	A.cols[i*dim*9 + 8*dim + d] = (i + nx + 1)*dim + d;
      }
    }
    else if (i==nx-1) {
      for (int d=0; d<2; d++){
	A.cols[i*dim*9 + 0*dim + d] = 0;
	A.cols[i*dim*9 + 1*dim + d] = 0;
	A.cols[i*dim*9 + 2*dim + d] = 0;
	A.cols[i*dim*9 + 3*dim + d] = (i - 1     )*dim + d;
	A.cols[i*dim*9 + 4*dim + d] = (i         )*dim + d;
	A.cols[i*dim*9 + 5*dim + d] = 0;
	A.cols[i*dim*9 + 6*dim + d] = (i + nx - 1)*dim + d;
	A.cols[i*dim*9 + 7*dim + d] = (i + nx    )*dim + d;
	A.cols[i*dim*9 + 8*dim + d] = 0;
      }
    }
    else if (i==nx*ny-1) {
      for (int d=0; d<2; d++){
	A.cols[i*dim*9 + 0*dim + d] = (i - nx - 1)*dim + d;
	A.cols[i*dim*9 + 1*dim + d] = (i - nx    )*dim + d;
	A.cols[i*dim*9 + 2*dim + d] = 0;
	A.cols[i*dim*9 + 3*dim + d] = (i - 1     )*dim + d;
	A.cols[i*dim*9 + 4*dim + d] = (i         )*dim + d;
	A.cols[i*dim*9 + 5*dim + d] = 0;
	A.cols[i*dim*9 + 6*dim + d] = 0;
	A.cols[i*dim*9 + 7*dim + d] = 0;
	A.cols[i*dim*9 + 8*dim + d] = 0;
      }
    }
    else if (i==(ny-1)*nx) {
      for (int d=0; d<2; d++){
	A.cols[i*dim*9 + 0*dim + d] = 0;
	A.cols[i*dim*9 + 1*dim + d] = (i - nx    )*dim + d;
	A.cols[i*dim*9 + 2*dim + d] = (i - nx + 1)*dim + d;
	A.cols[i*dim*9 + 3*dim + d] = 0;
	A.cols[i*dim*9 + 4*dim + d] = (i         )*dim + d;
	A.cols[i*dim*9 + 5*dim + d] = (i + 1     )*dim + d;
	A.cols[i*dim*9 + 6*dim + d] = 0;
	A.cols[i*dim*9 + 7*dim + d] = 0;
	A.cols[i*dim*9 + 8*dim + d] = 0;
      }
    }
  }
  cout << "A.num_rows = " << A.num_rows << endl;
  for (int i=0; i<(A.num_rows+1); i++)
    cout << A.row_offsets[i] << " ";
  cout << endl;
  for (int i=0; i<A.num_rows; i++) {
    for (int j=0; j<18; j++)
      cout << A.cols[i*18 + j] << " ";
    cout << endl;
  }

}

int csr_assembly_A (csr_matrix &A, Problem &problem)
{
  int dim = problem.dim;
  int npe = problem.npe;
  int nelem = problem.nelem; 
  int nx = problem.nx; 
  int ny = problem.ny; 

  int n0, n1, n2, n3;
  int cols_row_0[4]={4,5,7,8};
  int cols_row_1[4]={3,4,6,7};
  int cols_row_2[4]={0,1,3,4};
  int cols_row_3[4]={1,2,4,5};
  int yfactor;
  int xfactor;
  int index;

  double Ae[8][8];

  for (int e=0; e<nelem; e++){

    xfactor = e%nx;
    yfactor = e/ny;

    n0 = yfactor     * nx + xfactor     ;
    n1 = yfactor     * nx + xfactor + 1 ;
    n2 = (yfactor+1) * nx + xfactor + 1 ;
    n3 = (yfactor+1) * nx + xfactor     ;

    get_elemental (e, Ae);

      for (int n=0; n<4; n++){
	for (int d=0; d<2; d++){
	  A.coefs[n0*18 + cols_row_0[n]*2 + d    ] = Ae[0][n*2+d];
	  A.coefs[n0*18 + cols_row_0[n]*2 + d + 1] = Ae[1][n*2+d];
	  A.coefs[n1*18 + cols_row_1[n]*2 + d    ] = Ae[2][n*2+d];
	  A.coefs[n1*18 + cols_row_1[n]*2 + d + 1] = Ae[3][n*2+d];
	  A.coefs[n2*18 + cols_row_2[n]*2 + d    ] = Ae[4][n*2+d];
	  A.coefs[n2*18 + cols_row_2[n]*2 + d + 1] = Ae[5][n*2+d];
	  A.coefs[n3*18 + cols_row_3[n]*2 + d    ] = Ae[6][n*2+d];
	  A.coefs[n3*18 + cols_row_3[n]*2 + d + 1] = Ae[7][n*2+d];
	}
      }

  }

  return 0;
}

int main (int argc, char *argv[])
{
  Problem problem (argc, argv);
  csr_matrix A;
  csr_vector x, dx, res;

  double eps [] = { 0.005, 0.0, 0.0 };
  int nn = problem.nn;
  int dim = problem.dim;

  csr_alloc_A (A, nn*dim);
  csr_alloc_v (res, nn);
  csr_alloc_v (dx, nn);
  csr_alloc_v (x, nn);
  csr_set_A (A, problem);

  cout << "Init micropp -> preparing ..." << endl;
  cout << "dim = " << problem.dim << endl;
  cout << "nx  = " << problem.nx << endl;
  cout << "ny  = " << problem.ny << endl;
  cout << "nn  = " << problem.nn << endl;

  // assembly res

  // assembly A
  //  csr_assembly_A (A, problem);

  // solve system  A * dx = -res

  // free memory
  csr_free_A (A);
  csr_free_v (res);
  csr_free_v (dx);
  csr_free_v (x);

  return 0;
}

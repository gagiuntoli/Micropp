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

  A.row_offsets[0] = 0;
  for (int i=1; i<(A.num_rows+1); i++)
    A.row_offsets[i] =  18*i;
  
  for (int i=0; i<A.num_rows; i++){
    if (i>=nx && i<(ny-1)*nx) { // over y=0 and below y=ly
      if ((i%nx)>0 && (i%nx)<(nx-1)) {
	for (int d=0; d<2; d++){
	  A.cols[i*18 + 0*2 + d] = (i - nx - 1)*2 + d;
	  A.cols[i*18 + 1*2 + d] = (i - nx    )*2 + d;
	  A.cols[i*18 + 2*2 + d] = (i - nx + 1)*2 + d;
	  A.cols[i*18 + 3*2 + d] = (i - 1     )*2 + d;
	  A.cols[i*18 + 4*2 + d] = (i         )*2 + d;
	  A.cols[i*18 + 5*2 + d] = (i + 1     )*2 + d;
	  A.cols[i*18 + 6*2 + d] = (i + nx - 1)*2 + d;
	  A.cols[i*18 + 7*2 + d] = (i + nx    )*2 + d;
	  A.cols[i*18 + 8*2 + d] = (i + nx + 1)*2 + d;
	}
      }
    }
    else if (i==0) {       // the coorners
    }
    else if (i==nx-1) {
    }
    else if (i==nx*ny-1) {
    }
    else if (i==(ny-1)*nx) {
    }
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
  int n = 10;

  csr_alloc_A (A, n);
  csr_alloc_v (res, n);
  csr_alloc_v (dx, n);
  csr_alloc_v (x, n);
  csr_set_A (A, problem);

  cout << "Init micropp -> preparing ..." << endl;
  cout << "dim = " << problem.dim << endl;
  cout << "nx  = " << problem.nx << endl;
  cout << "ny  = " << problem.ny << endl;

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

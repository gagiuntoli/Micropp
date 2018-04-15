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
  cout << "A.num_rows = " << A.num_rows << endl;
  for (int i=0; i<(A.num_rows+1); i++)
    cout << A.row_offsets[i] << " ";
  cout << endl << endl;
  for (int i=0; i<A.num_rows; i++) {
    for (int j=0; j<18; j++)
      cout << A.cols[i*18 + j] << " ";
    cout << endl;
  }
  cout << endl ;

}

void get_elemental (int e, double (&Ae)[8][8])
{
  for (int i=0; i<8; i++)
    for (int j=0; j<8; j++)
      Ae[i][j] = 1;
}

int csr_assembly_A (csr_matrix &A, Problem &problem)
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

  for (int i=0; i<A.num_rows*18; i++)
    A.coefs[i] = 0.0;

  for (int e=0; e<nelem; e++){

    xfactor = e%(nx-1);
    yfactor = e/(ny-1);

    n0 = yfactor     * nx + xfactor     ;
    n1 = yfactor     * nx + xfactor + 1 ;
    n2 = (yfactor+1) * nx + xfactor + 1 ;
    n3 = (yfactor+1) * nx + xfactor     ;
    cout << "e : n0="<<n0<<" n1="<<n1<<" n2="<<n2<<" n3="<<n3<<endl;

    get_elemental (e, Ae);

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

  }

  for (int i=0; i<A.num_rows; i++) {
    for (int j=0; j<18; j++)
      cout << A.coefs[i*18 + j] << " ";
    cout << endl;
  }
  cout << endl ;

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
  cout << "dim   = " << problem.dim << endl;
  cout << "nx    = " << problem.nx << endl;
  cout << "ny    = " << problem.ny << endl;
  cout << "nn    = " << problem.nn << endl;
  cout << "nelem = " << problem.nelem << endl;

  // assembly res

  // assembly A
  csr_assembly_A (A, problem);

  // solve system  A * dx = -res

  // free memory
  csr_free_A (A);
  csr_free_v (res);
  csr_free_v (dx);
  csr_free_v (x);

  return 0;
}

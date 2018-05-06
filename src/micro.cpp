#include <vector>
#include <iostream>
#include "micro.h"

Problem::Problem (int dim, int size[3], int cg_its, double cg_tol)
{

  solver.max_its = cg_its;
  solver.min_tol = cg_tol;
  this->dim = dim;
  if (dim == 2) {
    nx = size[0];
    ny = size[1];
    nz = 1;
    nvoi = 3;
  }

  lx = 1.0;
  ly = 1.0;
  lz = 1.0;
  dx = lx/(nx-1);
  dy = ly/(ny-1);
  nn = nx*ny*nz;
  if (dim == 2) {
    nelem = (nx-1) * (ny-1);
  }

  ell_init (&A, nn*dim, nn*dim, 18);
  b  = (double*)malloc(nn*dim*sizeof(double));
  du = (double*)malloc(nn*dim*sizeof(double));
  u  = (double*)malloc(nn*dim*sizeof(double));
  stress = (double*)malloc(nelem*nvoi*sizeof(double));
  strain = (double*)malloc(nelem*nvoi*sizeof(double));
  elem_type = (int*)malloc(nelem*sizeof(int));
//  int_vars  = (double*)malloc(nn*dim*sizeof(double));

  for (int i=0; i<nn*dim; i++)
    u[i] = 0.0;

  for (int e=0; e<nelem; e++)
    elem_type[e] = getElemType(e);

  flag_print_A = false;
  flag_print_b = false;
  flag_print_u = false;
  flag_print_du = false;
  flag_print_newton = false;
  flag_print_solver = false;

  NewRap_Its = 3;
  NewRap_Tol = 1.0e-5;
}

int Problem::getElemType (int e)
{
  if (distance(e) < 0.2) {
    return 1;
  }
  else {
    return 0;
  }
  return -1;

}

//int csr_set_A (csr_matrix &A, Problem &problem)
//{
//  // set offsets and cols
//  int dim = problem.dim;
//  int npe = problem.npe;
//  int nelem = problem.nelem; 
//  int nx = problem.nx; 
//  int ny = problem.ny; 
//  int nn = problem.nn; 
//
//  A.row_offsets[0] = 0;
//  for (int i=1; i<(A.num_rows+1); i++)
//    A.row_offsets[i] =  18*i;
//  
//  for (int i=0; i<nn; i++){
//    // the coorners
//    if (i == 0) {       
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
//	}
//      }
//    }
//    else if (i == nx-1) {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
//	}
//      }
//    }
//    else if (i == nx*ny-1) {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
//	}
//      }
//    }
//    else if (i == (ny-1)*nx) {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
//	}
//      }
//    }
//    // y=0
//    else if (i < nx) {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
//	}
//      }
//    }
//    // y=ly
//    else if (i > (ny-1)*nx) {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
//	}
//      }
//    }
//    // x=0
//    else if ((i%nx) == 0) {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
//	}
//      }
//    }
//    // x=ly
//    else if ((i+1)%nx == 0) {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = 0; 
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = 0;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = 0;
//	}
//      }
//    }
//    // not y=0 neither y=ly, x=0, x=lx
//    else {
//      for (int d1=0; d1<2; d1++){
//	for (int d=0; d<2; d++){
//	  A.cols[i*dim*18 + 0*dim + d + 18*d1] = (i - nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 1*dim + d + 18*d1] = (i - nx    )*dim + d;
//	  A.cols[i*dim*18 + 2*dim + d + 18*d1] = (i - nx + 1)*dim + d;
//	  A.cols[i*dim*18 + 3*dim + d + 18*d1] = (i - 1     )*dim + d;
//	  A.cols[i*dim*18 + 4*dim + d + 18*d1] = (i         )*dim + d;
//	  A.cols[i*dim*18 + 5*dim + d + 18*d1] = (i + 1     )*dim + d;
//	  A.cols[i*dim*18 + 6*dim + d + 18*d1] = (i + nx - 1)*dim + d;
//	  A.cols[i*dim*18 + 7*dim + d + 18*d1] = (i + nx    )*dim + d;
//	  A.cols[i*dim*18 + 8*dim + d + 18*d1] = (i + nx + 1)*dim + d;
//	}
//      }
//    }
//  }
//
////  cout << "A.num_rows = " << A.num_rows << endl;
////  for (int i=0; i<(A.num_rows+1); i++)
////    cout << A.row_offsets[i] << " ";
////  cout << endl << endl;
////  for (int i=0; i<A.num_rows; i++) {
////    for (int j=0; j<18; j++)
////      cout << A.cols[i*18 + j] << " ";
////    cout << endl;
////  }
////  cout << endl ;
//
//}

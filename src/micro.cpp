#include <vector>
#include <boost/program_options.hpp>
#include <iostream>
#include "micro.h"

using namespace boost::program_options;

Problem::Problem (int argc, char *argv[])
{
  try {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("nx", value<int>()->default_value(10), "Num of Nodes in X dir")
      ("nr", value<int>()->default_value(10), "Num of Nodes in X dir")
      ("ny", value<int>()->default_value(10), "Num of Nodes in Y dir")
      ("nz", value<int>()->default_value(1) , "Num of Nodes in Z dir")
      ("dim", value<int>()->default_value(2), "Dimension");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help"))
      std::cout << desc << '\n';

    nx = vm["nx"].as<int>();
    ny = vm["ny"].as<int>();
    nz = vm["nz"].as<int>();
    dim = vm["dim"].as<int>();

  } catch (const error &ex) {
    std::cerr << ex.what() << '\n';
    throw 1;
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

  wg[0] = 0.25*(dx*dy);
  wg[1] = 0.25*(dx*dy);
  wg[2] = 0.25*(dx*dy);
  wg[3] = 0.25*(dx*dy);
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
/*   // for debugging
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

  double xg[4][2] = {
    {-0.577350269189626, -0.577350269189626},
    {+0.577350269189626, -0.577350269189626},
    {+0.577350269189626, +0.577350269189626},
    {-0.577350269189626, +0.577350269189626}};

  double dsh[4][2], b_mat[3][8], cxb[3][8];

  for (int i=0; i<8; i++) 
    for (int j=0; j<8; j++) 
      Ae[i][j] = 0.0;

  int dim = problem.dim;
  double dx = problem.dx;
  double dy = problem.dy;

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
	    Ae[i][j] += b_mat[m][i] * cxb[m][j] * problem.wg[gp];
	}
      }
    }
  } // gp loop
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

  for (int i=0; i<A.num_rows*18; i++)
    A.coefs[i] = 0.0;

  for (int e=0; e<nelem; e++) {

    xfactor = e%(nx-1);
    yfactor = e/(ny-1);

    n0 = yfactor     * nx + xfactor     ;
    n1 = yfactor     * nx + xfactor + 1 ;
    n2 = (yfactor+1) * nx + xfactor + 1 ;
    n3 = (yfactor+1) * nx + xfactor     ;
    //cout << "e : n0="<<n0<<" n1="<<n1<<" n2="<<n2<<" n3="<<n3<<endl;

    get_elemental (e, Ae, problem);

    for (int n=0; n<4; n++) {
      for (int d=0; d<2; d++) {
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

  /* boundary conditions */
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

#ifndef MICRO_HPP
#define MICRO_HPP

#include <vector>

class Problem {

  public:

    int npe = 4;
    int dim = 2;
    int nvoi= 3;

    std::vector<int> bc_nods, bc_y0, bc_y1, bc_x0, bc_x1;
    double xg[4][2];
    double wg[4];
    double b_mat[3][8][4];
    double dsh[4][2][4];
    int nx, ny, nz, nn; 
    double lx, ly, lz, dx, dy, dz;
    int nelem;
    int X0Y0_nod;
    int X1Y0_nod;
    int X1Y1_nod;
    int X0Y1_nod;
    int size_tot;
    int elem_type;
    std::vector<double> strain;
    std::vector<double> stress;
    std::vector<std::vector<int> > elements;
    std::vector<double> int_vars;
    double Ef, Em;
    double Sy_f, Sy_m;

    Problem (int argc, char *argv[]);

};

Problem::Problem (int argc, char *argv[])
{
  nx = 4;
  ny = 4;
  nz = 1;
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

  double xg[4][2] = {
    {-0.577350269189626, -0.577350269189626},
    {+0.577350269189626, -0.577350269189626},
    {+0.577350269189626, +0.577350269189626},
    {-0.577350269189626, +0.577350269189626}};

  for (int gp=0; gp<4; gp++) {
    dsh[0][0][gp] = -(1-xg[gp][1])/4*2/dx ; dsh[0][1][gp] = -(1-xg[gp][0])/4*2/dy;
    dsh[1][0][gp] = +(1-xg[gp][1])/4*2/dx ; dsh[1][1][gp] = -(1+xg[gp][0])/4*2/dy;
    dsh[2][0][gp] = +(1+xg[gp][1])/4*2/dx ; dsh[2][1][gp] = +(1+xg[gp][0])/4*2/dy;
    dsh[3][0][gp] = -(1+xg[gp][1])/4*2/dx ; dsh[3][1][gp] = +(1-xg[gp][0])/4*2/dy;
  }

  for (int gp=0; gp<4; gp++) {
    for (int i=0; i<4; i++) {
      b_mat[1][i*dim][gp] = dsh[i][0][gp]; b_mat[1][i*dim+1][gp] = 0            ;
      b_mat[2][i*dim][gp] = 0            ; b_mat[2][i*dim+1][gp] = dsh[i][1][gp];
      b_mat[3][i*dim][gp] = dsh[i][1][gp]; b_mat[3][i*dim+1][gp] = dsh[i][0][gp];
    }
  }
}

#endif

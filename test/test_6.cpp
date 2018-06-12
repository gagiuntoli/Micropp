#include <iostream>
#include <iomanip>
#include <ctime>
#include "micro.h"

using namespace std;

int main (int argc, char *argv[])
{

  int dim = 3;
  int nx = 10;
  int ny = 10;
  int nz = 10;
  double eps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  int size[3];
  size[0] = nx;
  size[1] = ny;
  size[2] = nz;

  int micro_type = 1; // 2 materiales matriz y fibra (3D esfera en matriz)
  double micro_params[4]; 
  micro_params[0] = 1.0; // lx
  micro_params[1] = 1.0; // ly
  micro_params[2] = 1.0; // lz
  micro_params[3] = 0.1; // grosor capa de abajo

  int mat_types[2]; // dos materiales lineales (type = 0)
  mat_types[0] = 1;
  mat_types[1] = 0;

  double mat_params[2*MAX_MAT_PARAM];
  mat_params[0*MAX_MAT_PARAM + 0] = 1.0e6; // E
  mat_params[0*MAX_MAT_PARAM + 1] = 0.3;   // nu
  mat_params[0*MAX_MAT_PARAM + 2] = 5.0e4; // Sy
  mat_params[0*MAX_MAT_PARAM + 3] = 5.0e4; // Ka

  mat_params[1*MAX_MAT_PARAM + 0] = 1.0e6;
  mat_params[1*MAX_MAT_PARAM + 1] = 0.3;
  mat_params[1*MAX_MAT_PARAM + 2] = 1.0e4;
  mat_params[1*MAX_MAT_PARAM + 3] = 0.0e-1;

  Problem micro (dim, size, micro_type, micro_params, mat_types, mat_params);

  int time_steps = 80;
  double stress_ave[6], ctan_ave[36];
  double d_eps = 0.01;
  int strain_comp = 2;

  for (int t=0; t<time_steps; t++) {

    cout << "Time step = " << t << endl;

    if (t<30)
      eps[strain_comp] += d_eps;
    else if (t<80)
      eps[strain_comp] -= d_eps;
    else if (t<130)
      eps[strain_comp] += d_eps;
    else if (t<250)
      eps[strain_comp] -= d_eps;
    else
      eps[strain_comp] += d_eps;

    double NR_norm;
    int NR_its, NR_non_linear;
    int LinCriteria;
    micro.loc_hom_Stress (1, eps, stress_ave);
    micro.getParams_NR (&NR_its, &NR_norm, &NR_non_linear);
    micro.getParams_LinCriteria (&LinCriteria);
    cout << "NEWTON-R ITS = " << NR_its << " TOL = " << NR_norm << " NON_LINEAR = " << NR_non_linear  << endl;
    cout << "LinCriteria = " << LinCriteria << endl;

    cout << "e11 = " << eps[strain_comp] << endl;
    cout 
      << "Average stress = " 
      << stress_ave[0] << " " << stress_ave[1] << " " << stress_ave[2] << " " 
      << stress_ave[3] << " " << stress_ave[4] << " " << stress_ave[5] 
      << endl;

    cout << endl;
//    micro.loc_hom_Ctan (1, eps, ctan_ave);
//
//    cout <<"Average Ctan = " << endl;
//    for (int i=0; i<6; i++) {
//      for (int j=0; j<6; j++)
//	cout << setw (8) << std::setprecision(3) << ctan_ave[i*6 + j] << " ";
//      cout << endl;
//    }

    micro.output (t, 1, 1, eps);
  }

  return 0;
}

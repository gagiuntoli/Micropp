#include <iostream>
#include <ctime>
#include <boost/program_options.hpp>
#include "micro.h"

using namespace std;
using namespace boost::program_options;

int main (int argc, char *argv[])
{

  try {

    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("nx", value<int>()->default_value(10), "Num of Nodes in X dir")
      ("ny", value<int>()->default_value(10), "Num of Nodes in Y dir")
      ("nz", value<int>()->default_value(10), "Num of Nodes in Z dir")
      ("dim", value<int>()->default_value(3), "Dimension");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help"))
      std::cout << desc << '\n';

    int dim = vm["dim"].as<int>();
    int nx = vm["nx"].as<int>();
    int ny = vm["ny"].as<int>();
    int nz = vm["nz"].as<int>();
    double eps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    int size[3];
    size[0] = nx;
    size[1] = ny;
    size[2] = nz;

    int micro_type = 0; // 2 materiales matriz y fibra (3D esfera en matriz)
    double micro_params[4]; 
    micro_params[0] = 1.0; // lx
    micro_params[1] = 1.0; // ly
    micro_params[2] = 1.0; // lz
    micro_params[3] = 0.2; // radio de la esfera

    int mat_types[2]; // dos materiales lineales (type = 0)
    mat_types[0] = 1;
    mat_types[1] = 0;

    double params[2*MAX_MAT_PARAM];
    params[0*MAX_MAT_PARAM + 0] = 1.0e6;
    params[0*MAX_MAT_PARAM + 1] = 0.3;
    params[0*MAX_MAT_PARAM + 2] = 1.0e4;
    params[0*MAX_MAT_PARAM + 3] = 1.0e4;

    params[1*MAX_MAT_PARAM + 0] = 1.0e7;
    params[1*MAX_MAT_PARAM + 1] = 0.3;
    params[1*MAX_MAT_PARAM + 2] = 1.0e4;
    params[1*MAX_MAT_PARAM + 3] = 1.0e2;

    Problem micro (dim, size, micro_type, micro_params, mat_types, params);

    int time_steps = 10;
    double stress_ave[6];

    for (int t=0; t<time_steps; t++) {

      cout << "Time step = " << t << endl;

      eps[0] = t*1.0/time_steps * 0.002;
      micro.loc_hom_Stress (1, eps, stress_ave);

      cout 
	<< "Average stress = " 
	<< stress_ave[0] << " " << stress_ave[1] << " " << stress_ave[2] << " " 
	<< stress_ave[3] << " " << stress_ave[4] << " " << stress_ave[5] 
	<< endl;

      micro.output (t, 1, 1, eps);
    }

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}

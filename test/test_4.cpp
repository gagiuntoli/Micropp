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
      ("cg_its", value<int>()->default_value(1000), "Max Number of iterations (CG)")
      ("cg_tol", value<double>()->default_value(1e-8), "Minimun Tolerance (CG)")
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
    int cg_its = vm["cg_its"].as<int>();
    double cg_tol = vm["cg_tol"].as<double>();
    double eps[6] = {0.005, 0.0, 0.0, 0.0, 0.0, 0.0};

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

    int types[2]; // dos materiales lineales (type = 0)
    types[0] = 0;
    types[1] = 0;

    double params[2*MAX_MAT_PARAM];
    params[0*MAX_MAT_PARAM + 0] = 1.0e6;
    params[0*MAX_MAT_PARAM + 1] = 0.3;

    params[1*MAX_MAT_PARAM + 0] = 1.0e7;
    params[1*MAX_MAT_PARAM + 1] = 0.3;

    Problem micro (dim, size, micro_type, micro_params, types, params);

    double start, end, t_lap;

    start = clock();
    micro.setDisp(eps);
    bool non_linear;
    double norm = micro.Assembly_b(NULL, NULL, &non_linear);
    cout << "Residual norm = " << norm << endl;
    micro.Assembly_A(NULL);
    end = clock();
    t_lap = double(end - start) / CLOCKS_PER_SEC;
    cout << "time assembly : " << t_lap << endl;

    start = clock();
    micro.solve();
    end = clock();
    t_lap = double(end - start) / CLOCKS_PER_SEC;
    cout << "time solve : " << t_lap << endl;

    for (int i=0; i<micro.nn*dim; i++)
      micro.u[i] = micro.u[i] + micro.du[i];

    norm = micro.Assembly_b(NULL, NULL, &non_linear);
    cout << "Residual norm = " << norm << endl;

    micro.calcDistributions(NULL);
    micro.writeVtu (1, 2, NULL);

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}

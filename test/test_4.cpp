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
    double eps[3] = {0.005, 0.0, 0.0};

    int size[3];
    size[0] = nx;
    size[1] = ny;
    size[2] = nz;

    double start, end, t_lap;
    Problem micro (dim, size, cg_its, cg_tol);

    micro.setDisp(eps);

    micro.writeVtu (1, 2);

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}

#include <iostream>
#include <ctime>
#include <boost/program_options.hpp>
#include "micro.h"

using namespace std;
using namespace boost::program_options;

int main (int argc, char *argv[])
{

  try {

    bool flag_print_A=false;
    bool flag_print_b=false;
    bool flag_print_u=false;
    bool flag_print_du=false;

    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("nx", value<int>()->default_value(10), "Num of Nodes in X dir")
      ("nr", value<int>()->default_value(10), "Num of Nodes in X dir")
      ("ny", value<int>()->default_value(10), "Num of Nodes in Y dir")
      ("nz", value<int>()->default_value(1) , "Num of Nodes in Z dir")
      ("cg_its", value<int>()->default_value(100)     , "Max Number of iterations (CG)")
      ("cg_tol", value<double>()->default_value(1e-8) , "Minimun Tolerance (CG)")
      ("print_A", bool_switch(&flag_print_A), "Print ELL Matrix")
      ("print_b", bool_switch(&flag_print_b), "Print RHS b")
      ("print_u", bool_switch(&flag_print_u), "Print Solution u")
      ("print_du", bool_switch(&flag_print_du), "Print Increment du")
      ("dim", value<int>()->default_value(2), "Dimension");

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

    int size[3];
    size[0] = nx;
    size[1] = ny;
    size[2] = nz;

    double start, end, t_lap;
    Problem micro (dim, size, cg_its, cg_tol);

    micro.flag_print_A = flag_print_A;
    micro.flag_print_b = flag_print_b;
    micro.flag_print_u = flag_print_u;
    micro.flag_print_du = flag_print_du;

    // assembly
    start = clock();
    double eps[3] = {0.005, 0.0, 0.0};
    micro.setDisp(eps);
    micro.Assembly_A();
    micro.Assembly_b();
    end = clock();
    t_lap = double(end - start) / CLOCKS_PER_SEC;
    cerr << "Time Assembly : " << t_lap << endl;

    // solve
    start = clock();
    micro.solve();
    end = clock();
    t_lap = double(end - start) / CLOCKS_PER_SEC;
    cerr << "Time Solving  : " << t_lap << endl;

    // writting
    start = clock();
    micro.write_vtu();
    end = clock();
    t_lap = double(end - start) / CLOCKS_PER_SEC;
    cerr << "Time Writing  : " << t_lap << endl;

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}

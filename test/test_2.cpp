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
      ("nr", value<int>()->default_value(10), "Num of Nodes in X dir")
      ("ny", value<int>()->default_value(10), "Num of Nodes in Y dir")
      ("nz", value<int>()->default_value(1) , "Num of Nodes in Z dir")
      ("cg_its", value<int>()->default_value(1000)    , "Max Number of iterations (CG)")
      ("cg_tol", value<double>()->default_value(1e-8) , "Minimun Tolerance (CG)")
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
    double eps[3] = {0.005, 0.0, 0.0};

    int size[3];
    size[0] = nx;
    size[1] = ny;
    size[2] = nz;

    double start, end, t_lap;
    Problem micro (dim, size, cg_its, cg_tol);

    // test localization and homogenization of stress
    double stress_mac[6];
    start = clock();
    micro.loc_hom_Stress(eps, stress_mac);
    cout << "The average stress for loc-hom is : " << stress_mac[0] << " " << stress_mac[1] << " " << stress_mac[2] << endl;
    end = clock();
    t_lap = double(end - start) / CLOCKS_PER_SEC;
    cout << "Time calc MacroStress: " << t_lap << endl;

    // test MacroCtan calculation
    double MacroCtan[81];
    start = clock();
    micro.loc_hom_Ctan (eps, MacroCtan);
    cout << "MacroCtan: " << endl;
    for (int i=0; i<micro.nvoi; i++){
      for (int j=0; j<micro.nvoi; j++)
	cout << MacroCtan[i*micro.nvoi+j] << " ";
      cout << endl;
    }
    end = clock();
    t_lap = double(end - start) / CLOCKS_PER_SEC;
    cout << "Time calc MacroCtan: " << t_lap << endl;

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}

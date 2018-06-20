/*
 *  MicroPP : finite element library to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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
    double eps1[6] = {0.005, 0.0, 0.0, 0.0, 0.0, 0.0};

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
    mat_types[0] = 0;
    mat_types[1] = 0;

    double params[2*MAX_MAT_PARAM];
    params[0*MAX_MAT_PARAM + 0] = 1.0e6;
    params[0*MAX_MAT_PARAM + 1] = 0.3;
    params[0*MAX_MAT_PARAM + 2] = 1.0e4;
    params[0*MAX_MAT_PARAM + 3] = 1.0e2;

    params[1*MAX_MAT_PARAM + 0] = 1.0e7;
    params[1*MAX_MAT_PARAM + 1] = 0.3;
    params[1*MAX_MAT_PARAM + 2] = 1.0e4;
    params[1*MAX_MAT_PARAM + 3] = 1.0e2;

    Problem micro1 (dim, size, micro_type, micro_params, mat_types, params);

    bool non_linear;
    micro1.setDisp(eps1);
    micro1.newtonRaphson(&non_linear);

    double stress_ave[6];
    micro1.calcAverageStress(stress_ave);
    cout 
      << "Average stress = " 
      << stress_ave[0] << " " << stress_ave[1] << " " << stress_ave[2] << " " 
      << stress_ave[3] << " " << stress_ave[4] << " " << stress_ave[5] 
      << endl;

    double strain_ave[6];
    micro1.calcAverageStrain(strain_ave);
    cout 
      << "Average strain = " 
      << strain_ave[0] << " " << strain_ave[1] << " " << strain_ave[2] << " " 
      << strain_ave[3] << " " << strain_ave[4] << " " << strain_ave[5]
      << endl;

    micro1.calcDistributions();
    micro1.writeVtu (1, 2);

    double eps2[6] = {0.0, 0.0, 0.0, 0.005, 0.0, 0.0};
    micro_type = 1; // 2 materiales matriz y fibra (3D 2 capas planas)
    micro_params[0] = 1.0; // lx
    micro_params[1] = 1.0; // ly
    micro_params[2] = 1.0; // lz
    micro_params[3] = 0.2; // espesor de la capa 0

    mat_types[0] = 0;
    mat_types[1] = 0;

    params[0*MAX_MAT_PARAM + 0] = 1.0e6;
    params[0*MAX_MAT_PARAM + 1] = 0.3;

    params[1*MAX_MAT_PARAM + 0] = 1.0e7;
    params[1*MAX_MAT_PARAM + 1] = 0.3;

    Problem micro2 (dim, size, micro_type, micro_params, mat_types, params);

    micro2.setDisp(eps2);
    micro2.newtonRaphson(&non_linear);

    micro2.calcAverageStress(stress_ave);
    cout 
      << "Average stress = " 
      << stress_ave[0] << " " << stress_ave[1] << " " << stress_ave[2] << " " 
      << stress_ave[3] << " " << stress_ave[4] << " " << stress_ave[5] 
      << endl;

    micro2.calcAverageStrain(strain_ave);
    cout 
      << "Average strain = " 
      << strain_ave[0] << " " << strain_ave[1] << " " << strain_ave[2] << " " 
      << strain_ave[3] << " " << strain_ave[4] << " " << strain_ave[5]
      << endl;

    micro2.calcDistributions();
    micro2.writeVtu (2, 2);

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}

/*
 *  This is a test example for MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
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
#include <iomanip>

#include <ctime>
#include <cassert>

#include "material.hpp"

using namespace std;


void micropp_simulator(material_t *materials)
{
	cout << endl << "MicroPP Simulator" << endl;

	material_t *micropp_material_list[2];

	for (int i = 0; i < 2; ++i) {
		double params[3];
		micropp_material_list[i] = material_t::make_material(materials[i]);
	}
	for (auto it : micropp_material_list)
		it->print_n();

	double eps[] = { 0.2, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double sig[] = { 0.2, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double *vars = nullptr;
	for (auto it : micropp_material_list) {
		it->get_stress(eps, sig, vars);
		for (int i = 0; i < 6; ++i)
			cout << " " << sig[i];
		cout << endl;
	}
}

int main (void)
{

	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.3, 5.0e4, 2.0e4, 0);
	mat_params[1].set(1.0e3, 0.3, 5.0e4, 1.0e3, 2);

	micropp_simulator(mat_params);

	material_t *mat_params_n[2];
	mat_params_n[0] = new material_elastic(1.0, 2.0);
	mat_params_n[1] = new material_damage(1.0, 2.0, 3.0);

	for (auto it : mat_params_n)
		it->print_n();

	return 0;
}

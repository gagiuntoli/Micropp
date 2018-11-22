/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                        Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef MICROPP_C_WRAPPER_H
#define MICROPP_C_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

	void micropp_C_material_set(int num_mat, const double E, double nu,
				    double Ka, double Sy, int type);
	void micropp_C_material_print(int num_mat);

	void micropp_C_create3(int ngp, int size[3], int type, double *params);
	void micropp_C_destroy3();
	void micropp_C_set_strain3(int gp, double strain[6]);
	void micropp_C_get_stress3(int gp, double stress[6]);
	void micropp_C_get_ctan3(int gp, double ctan[36]);
	int micropp_C_get_non_linear_gps(void);
	void micropp_C_update_vars();
	void micropp_C_homogenize();
	void micropp_C_print_info();
	int micropp_C_get_sigma_cost3(int gp);

#ifdef __cplusplus
}
#endif
#endif

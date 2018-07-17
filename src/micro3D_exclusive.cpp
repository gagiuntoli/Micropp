/*
 *  This source code is part of MicroPP: a finite element library
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

#include <cmath>
#include <iostream>

#include "instrument.hpp"
#include "micro.hpp"

#define nod_index(i,j,k) ((k)*nx*ny + (j)*nx + (i))

using namespace std;

template <>
void micropp<3>::get_ctan_plast_sec(int ex, int ey, int ez,
                                    int gp, double ctan[6][6]) const
{
	INST_START;
	bool non_linear;
	double stress_pert[6], strain_pert[6], strain_0[6], d_strain = 1.0e-8;
	get_strain(gp, strain_0, ex, ey, ez);

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++)
			strain_pert[j] = 0.0;

		if (fabs(strain_0[i]) > 1.0e-7)
			strain_pert[i] = strain_0[i];
		else
			strain_pert[i] = d_strain;

		get_stress(gp, strain_pert, &non_linear, stress_pert, ex, ey, ez);

		for (int j = 0; j < 6; j++)
			ctan[j][i] = stress_pert[j] / strain_pert[i];
	}

}


template <>
void micropp<3>::get_ctan_plast_exact(int ex, int ey, int ez,
                                      int gp, double ctan[6][6]) const
{
	INST_START;
	double strain[6];
	int e = glo_elem3D(ex, ey, ez);
	get_strain(gp, strain, ex, ey, ez);

	const material_t material = get_material(e);

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			ctan[i][j] = 0.0;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ctan[i][j] += material.lambda;

	for (int i = 0; i < 3; i++)
		ctan[i][i] += 2 * material.mu;

	for (int i = 3; i < 6; i++)
		ctan[i][i] += material.mu;

}

template <>
void micropp<3>::get_ctan_plast_pert(int ex, int ey, int ez,
                                     int gp, double ctan[6][6]) const
{
	INST_START;
	bool non_linear;
	double stress_0[6], stress_pert[6], strain_0[6], strain_pert[6], deps = 1.0e-8;
	get_strain(gp, strain_0, ex, ey, ez);
	get_stress(gp, strain_0, &non_linear, stress_0, ex, ey, ez);

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++)
			strain_pert[j] = strain_0[j];

		strain_pert[i] += deps;
		get_stress(gp, strain_pert, &non_linear, stress_pert, ex, ey, ez);

		for (int j = 0; j < 6; j++)
			ctan[j][i] = (stress_pert[j] - stress_0[j]) / deps;
	}
}


template <>
void micropp<3>::get_dev_tensor(double tensor[6], double tensor_dev[6]) const
{
	for (int i = 0; i < 6; i++)
		tensor_dev[i] = tensor[i];
	for (int i = 0; i < 3; i++)
		tensor_dev[i] -= (1 / 3.0) * (tensor[0] + tensor[1] + tensor[2]);
}


template <>
void micropp<3>::plastic_step(const material_t *material, double eps[6],
                             double eps_p_old[6], double alpha_old,
                             double eps_p_new[6], double *alpha_new,
                             bool *non_linear, double stress[6]) const
{
	double eps_dev[6];
	double eps_p_dev_1[6];
	double normal[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double dl = 0.0;
	double sig_dev_trial[6], sig_dev_trial_norm;

	get_dev_tensor(eps_p_old, eps_p_dev_1);
	get_dev_tensor(eps, eps_dev);

	for (int i = 0; i < 3; ++i)
		sig_dev_trial[i] = 2 * material->mu * (eps_dev[i] - eps_p_dev_1[i]);

	for (int i = 3; i < 6; ++i)
		sig_dev_trial[i] = material->mu * (eps_dev[i] - eps_p_dev_1[i]);

	sig_dev_trial_norm = sqrt(sig_dev_trial[0] * sig_dev_trial[0] +
	                          sig_dev_trial[1] * sig_dev_trial[1] +
	                          sig_dev_trial[2] * sig_dev_trial[2] +
	                          2 * sig_dev_trial[3] * sig_dev_trial[3] +
	                          2 * sig_dev_trial[4] * sig_dev_trial[4] +
	                          2 * sig_dev_trial[5] * sig_dev_trial[5]);

	double f_trial = sig_dev_trial_norm - sqrt(2.0 / 3) * (material->Sy + material->Ka * alpha_old);

	if (f_trial > 0) {
		*non_linear = true;

		for (int i = 0; i < 6; ++i)
			normal[i] = sig_dev_trial[i] / sig_dev_trial_norm;

		dl = f_trial / (2 * material->mu * (1.0 + (0.0 * material->Ka) / (3 * material->mu)));

		if (eps_p_new && alpha_new) {
			for (int i = 0; i < 6; ++i)
				eps_p_new[i] = eps_p_old[i] + dl * normal[i];
			*alpha_new = alpha_old + sqrt(2.0 / 3) * dl;
		}
	} else if (eps_p_new && alpha_new){
		for (int i = 0; i < 6; ++i)
			eps_p_new[i] = eps_p_old[i];
		*alpha_new = alpha_old;
	}

	//sig_2 = s_trial + K * tr(eps) * 1 - 2*mu*dl*normal;
	for (int i = 0; i < 6; ++i)
		stress[i] = sig_dev_trial[i];

	for (int i = 0; i < 3; ++i)
		stress[i] += material->k * (eps[0] + eps[1] + eps[2]);

	for (int i = 0; i < 6; ++i)
		stress[i] -= 2 * material->mu * dl * normal[i];
}



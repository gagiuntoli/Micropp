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
 *  
 *************************************************************
 *
 * This code was though only for GPUs do to the current incompatibility
 * between OpenACC and the "virtual functions".
 *
 */


#include "material_acc.hpp"


material_acc::material_acc(const material_base material)
{
		E = material.E;
		nu = material.nu;
		Ka = material.Ka;
		Sy = material.Sy;
		k = material.k;
		mu = material.mu;
		lambda = material.lambda;
		Xt = material.Xt;
		type = material.type;
}


void material_acc::get_stress(const double *eps, double *stress,
					const double *history_params) const
{
	if(type == MATERIAL_ELASTIC) {
		get_stress_elastic(this, eps, stress, history_params);
	}
}


void material_acc::get_ctan(const double *eps, double *ctan,
				      const double *history_params) const
{
	if(type == MATERIAL_ELASTIC) {
		get_ctan_elastic(this, eps, ctan, history_params);
	}
}


void apply_perturbation(const material_acc *material,
			const double *eps, double *ctan,
			const double *vars_old,
			void (*get_stress_ptr)(const material_acc *, const double *, double *, const double *))
{
	double stress_0[6];
	(*get_stress_ptr)(material, eps, stress_0, vars_old);

	for (int i = 0; i < 6; ++i) {

		double eps_1[6];
		memcpy(eps_1, eps, 6 * sizeof(double));
		eps_1[i] += D_EPS_CTAN;

		double stress_1[6];
		(*get_stress_ptr)(material, eps_1, stress_1, vars_old);

		for (int j = 0; j < 6; ++j)
			ctan[j * 6 + i] =
				(stress_1[j] - stress_0[j]) / D_EPS_CTAN;
	}
}


static void get_dev_tensor(const double tensor[6], double tensor_dev[6])
{
	memcpy(tensor_dev, tensor, 6 * sizeof(double));
	for (int i = 0; i < 3; i++)
		tensor_dev[i] -= \
			(1 / 3.0) * (tensor[0] + tensor[1] + tensor[2]);
}


// ELASTIC MATERIAL


void get_stress_elastic(const material_acc *material,
			const double *eps, double *stress,
			const double *history_params)
{
	// stress[i][j] = lambda eps[k][k] * delta[i][j] + mu eps[i][j]
	for (int i = 0; i < 3; ++i)
		stress[i] = material->lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * material->mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = material->mu * eps[i];
}


void get_ctan_elastic(const material_acc *material,
		      const double *eps, double *ctan,
		      const double *history_params)
{
	// C = lambda * (1x1) + 2 mu I
	memset(ctan, 0, 6 * 6 * sizeof(double));

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			ctan[i * 6 + j] += material->lambda;

	for (int i = 0; i < 3; ++i)
		ctan[i * 6 + i] += 2 * material->mu;

	for (int i = 3; i < 6; ++i)
		ctan[i * 6 + i] = material->mu;
}


bool evolute_elastic(const double *eps, const double *vars_old,
		     double *vars_new)
{
	// we don't have to evolute nothing is always linear
	return false;
}


// PLASTIC MATERIAL


bool plastic_law(const material_acc *material, const double eps[6],
		 const double *_eps_p_old, const double *_alpha_old,
		 double *_dl, double _normal[6], double _s_trial[6])
{
	/*
	 * Calculates _dl, _normal and _s_trial to used in the other plastic 
	 * material functions. Returns <true> if it enters in non-linear zone
	 * <false> if not.
	 */

	const double zeros[6] = { 0.0 };
	const double alpha_old = (_alpha_old) ? *_alpha_old : 0;
	const double *eps_p_old = (_eps_p_old) ? _eps_p_old : zeros;

	double eps_dev[6], eps_p_dev_1[6];

	get_dev_tensor(eps_p_old, eps_p_dev_1);
	get_dev_tensor(eps, eps_dev);

	for (int i = 0; i < 3; ++i)
		_s_trial[i] = 2 * material->mu * (eps_dev[i] - eps_p_dev_1[i]);

	for (int i = 3; i < 6; ++i)
		_s_trial[i] = material->mu * (eps_dev[i] - eps_p_dev_1[i]);

	double tmp = 0.0;
	for (int i = 0; i < 6; ++i)
		tmp += _s_trial[i] * _s_trial[i];
	double s_norm = sqrt(tmp);

	double f_trial = s_norm - SQRT_2DIV3 * (material->Sy + material->Ka * alpha_old);

	if (f_trial > 0) {

		for (int i = 0; i < 6; ++i)
			_normal[i] = _s_trial[i] / s_norm;
		*_dl = f_trial / (2. * material->mu * (1. + material->Ka / (3. * material->mu)));
		return true;

	} else {

		memset(_normal, 0, 6 * sizeof(double));
		*_dl = 0;
		return false;
	}
}


void get_stress_plastic(const material_acc *material,
			const double *eps, double *stress,
			const double *history_params)
{
	double dl, normal[6], s_trial[6];
	const double *eps_p_old = (history_params) ? &(history_params[0]) : nullptr;
	const double *alpha_old = (history_params) ? &(history_params[6]) : nullptr;

	bool nl_flag = plastic_law(material, eps, eps_p_old, alpha_old, &dl, normal,
				   s_trial);

	//sig_2 = s_trial + K * tr(eps) * 1 - 2 * mu * dl * normal;
	memcpy(stress, s_trial, 6 * sizeof(double));

	for (int i = 0; i < 3; ++i)
		stress[i] += material->k * (eps[0] + eps[1] + eps[2]);

	for (int i = 0; i < 6; ++i)
		stress[i] -= 2 * material->mu * dl * normal[i];
}


void get_ctan_plastic(const material_acc *material, const double *eps,
		      double *ctan, const double *vars_old)
{
	apply_perturbation(material, eps, ctan, vars_old, get_stress_plastic);
}


bool evolute_plastic(const material_acc *material, const double *eps,
		     const double *vars_old, double *vars_new)
{
	const double *eps_p_old = (vars_old) ? &(vars_old[0]) : nullptr;
	const double *alpha_old = (vars_old) ? &(vars_old[6]) : nullptr;
	double *eps_p_new = (vars_new) ? &(vars_new[0]) : nullptr;
	double *alpha_new = (vars_new) ? &(vars_new[6]) : nullptr;

	double dl, normal[6], s_trial[6];
	bool nl_flag = plastic_law(material, eps, eps_p_old, alpha_old, &dl, normal,
				   s_trial);

	if (eps_p_old != nullptr && eps_p_new != nullptr)
		for (int i = 0; i < 6; ++i)
			eps_p_new[i] = eps_p_old[i] + dl * normal[i];

	if (alpha_old != nullptr && alpha_new != nullptr)
		*alpha_new =  *alpha_old + SQRT_2DIV3 * dl + 0;

	return nl_flag;
}


// DAMAGE MATERIAL


bool damage_law(const material_acc *material,
		const double *eps, const double e_old,
		const double D_old, double *_e, 
		double *_D, double *stress)
{
	/*
	 * Calculates the linear stree <stress>, and <e> and <D> using the
	 * strain <eps> and <e_old>
	 *
	 * e_old = vars_old[0]
	 *
	 */

	// First suppose we are in linear zone
	for (int i = 0; i < 3; ++i)
		stress[i] = material->lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * material->mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = material->mu * eps[i];

	// Now check if we have entered in non-linear zone
	double e = 0.0;
	for (int i = 0; i < 3; ++i)
		e += stress[i] * stress[i];
	e /= (material->Xt * material->Xt);
	e = sqrt(e);

	e = max(e_old, e);

	double D = (e < 1.0) ? 0.0 : (1 - exp(1 - e));
	D = min(D, D_old + 0.02);

	if (_e != nullptr)
		*_e = e;

	if (_D != nullptr)
		*_D = D;

	return ((e < 1.0) ? false : true);
}


void get_stress_damage(const material_acc *material,
		       const double *eps, double *stress,
		       const double *vars_old)
{
	/*
	 * Calculates the <stress> according to <eps> and <vars_old>.
	 *
	 * e_old = vars_old[0]
	 *
	 */

	const double e_old = (vars_old != nullptr) ? vars_old[0] : 0.0;
	const double D_old = (vars_old != nullptr) ? vars_old[1] : 0.0;
	double D;
	damage_law(material, eps, e_old, D_old, nullptr, &D, stress);

	for (int i = 0; i < 6; ++i)
		stress[i] *= (1 - D);
}


void get_ctan_damage(const material_acc *material,
		     const double *eps, double *ctan,
		     const double *vars_old)
{
	apply_perturbation(material, eps, ctan, vars_old, get_stress_damage);
}


bool evolute_damage(const material_acc *material,
		    const double *eps, const double *vars_old,
		    double *vars_new)
{
	/* Assign new values to <vars_new> according to <eps> and <vars_old>.
	 * returns <true> if the material has entered in non-linear range, 
	 * <false> if not.
	 */
	const double e_old = (vars_old) ? vars_old[0] : 0;
	const double D_old = (vars_old) ? vars_old[1] : 0;
	double *e_new = (vars_new) ? &(vars_new[0]) : nullptr;
	double *D_new = (vars_new) ? &(vars_new[1]) : nullptr;

	double stress[6];

	bool non_linear = damage_law(material, eps, e_old, D_old, e_new, D_new, stress);

	return non_linear;
}

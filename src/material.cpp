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


#include "material.hpp"


material_t *material_t::make_material(const struct material_base material)
{
	/*
	 * This fabric creates the corresponding subclass according to the 
	 * data members of <material>. This last is a "struct" than can be
	 * manage also from C/Fortran easily and that is why we passed micropp
	 * constructor the data in this form.
	 */

	switch (material.type) {
		case 0:
			return new material_elastic(material.E, material.nu);
			break;
		case 1:
			return new material_plastic(material.E, material.nu,
						    material.Ka, material.Sy);
			break;
		case 2:
			return new material_damage(material.E, material.nu,
						   material.Xt);
			break;
		default:
			break;
	}
}


void material_t::apply_perturbation(const double *eps, double *ctan,
				    const double *vars_old) const
{
	double stress_0[6];
	get_stress(eps, stress_0, vars_old);

	for (int i = 0; i < 6; ++i) {

		double eps_1[6];
		memcpy(eps_1, eps, 6 * sizeof(double));
		eps_1[i] += D_EPS_CTAN;

		double stress_1[6];
		get_stress(eps_1, stress_1, vars_old);

		for (int j = 0; j < 6; ++j)
			ctan[j * 6 + i] =
				(stress_1[j] - stress_0[j]) / D_EPS_CTAN;
	}
}


void get_dev_tensor(const double tensor[6], double tensor_dev[6])
{
	memcpy(tensor_dev, tensor, 6 * sizeof(double));
	for (int i = 0; i < 3; i++)
		tensor_dev[i] -= \
			(1 / 3.0) * (tensor[0] + tensor[1] + tensor[2]);
}


// ELASTIC MATERIAL


void material_elastic::get_stress(const double *eps, double *stress,
				  const double *history_params) const
{
	// stress[i][j] = lambda eps[k][k] * delta[i][j] + mu eps[i][j]
	for (int i = 0; i < 3; ++i)
		stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = mu * eps[i];
}


void material_elastic::get_ctan(const double *eps, double *ctan,
				const double *history_params) const
{
	// C = lambda * (1x1) + 2 mu I
	memset(ctan, 0, 6 * 6 * sizeof(double));

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			ctan[i * 6 + j] += lambda;

	for (int i = 0; i < 3; ++i)
		ctan[i * 6 + i] += 2 * mu;

	for (int i = 3; i < 6; ++i)
		ctan[i * 6 + i] = mu;
}


bool material_elastic::evolute(const double *eps, const double *vars_old,
			       double *vars_new) const
{
	// we don't have to evolute nothing is always linear
	return false;
}


void material_elastic::print() const
{
	cout << "Type : Elastic" << endl;
	cout << scientific << "E = " << E << " nu = " << nu << endl;
}


// PLASTIC MATERIAL


bool material_plastic::plastic_law(const double eps[6],
				   const double *_eps_p_old,
				   const double *_alpha_old,
				   double *_dl,
				   double _normal[6],
				   double _s_trial[6]) const
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
		_s_trial[i] = 2 * mu * (eps_dev[i] - eps_p_dev_1[i]);

	for (int i = 3; i < 6; ++i)
		_s_trial[i] = mu * (eps_dev[i] - eps_p_dev_1[i]);

	double tmp = 0.0;
	for (int i = 0; i < 6; ++i)
		tmp += _s_trial[i] * _s_trial[i];
	double s_norm = sqrt(tmp);

	double f_trial = s_norm - SQRT_2DIV3 * (Sy + Ka * alpha_old);

	if (f_trial > 0) {

		for (int i = 0; i < 6; ++i)
			_normal[i] = _s_trial[i] / s_norm;
		*_dl = f_trial / (2. * mu * (1. + Ka / (3. * mu)));
		return true;

	} else {

		memset(_normal, 0, 6 * sizeof(double));
		*_dl = 0;
		return false;
	}
}


void material_plastic::get_stress(const double *eps, double *stress,
				  const double *history_params) const
{
	double dl, normal[6], s_trial[6];
	const double *eps_p_old = (history_params) ? &(history_params[0]) : nullptr;
	const double *alpha_old = (history_params) ? &(history_params[6]) : nullptr;

	bool nl_flag = plastic_law(eps, eps_p_old, alpha_old, &dl, normal,
				   s_trial);

	//sig_2 = s_trial + K * tr(eps) * 1 - 2 * mu * dl * normal;
	memcpy(stress, s_trial, 6 * sizeof(double));

	for (int i = 0; i < 3; ++i)
		stress[i] += k * (eps[0] + eps[1] + eps[2]);

	for (int i = 0; i < 6; ++i)
		stress[i] -= 2 * mu * dl * normal[i];
}


void material_plastic::get_ctan(const double *eps, double *ctan,
				const double *vars_old) const
{
	apply_perturbation(eps, ctan, vars_old);
}


bool material_plastic::evolute(const double *eps, const double *vars_old,
			       double *vars_new) const
{
	const double *eps_p_old = (vars_old) ? &(vars_old[0]) : nullptr;
	const double *alpha_old = (vars_old) ? &(vars_old[6]) : nullptr;
	double *eps_p_new = (vars_new) ? &(vars_new[0]) : nullptr;
	double *alpha_new = (vars_new) ? &(vars_new[6]) : nullptr;

	double dl, normal[6], s_trial[6];
	bool nl_flag = plastic_law(eps, eps_p_old, alpha_old, &dl, normal,
				   s_trial);

	if (eps_p_old != nullptr && eps_p_new != nullptr)
		for (int i = 0; i < 6; ++i)
			eps_p_new[i] = eps_p_old[i] + dl * normal[i];

	if (alpha_old != nullptr && alpha_new != nullptr)
		*alpha_new =  *alpha_old + SQRT_2DIV3 * dl + 0;

	return nl_flag;
}


void material_plastic::print() const
{
	cout << "Type : Plastic" << endl;
	cout << "E = " << E << " nu = " << nu << " Ka = " << Ka << " Sy = "
		<< Sy << endl;
}


// DAMAGE MATERIAL


bool material_damage::damage_law(const double *eps, const double e_old,
				 const double D_old, double *_e, 
				 double *_D, double *stress) const
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
		stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = mu * eps[i];

	// Now check if we have entered in non-linear zone
	double e = 0.0;
	for (int i = 0; i < 3; ++i)
		e += stress[i] * stress[i];
	e /= (Xt * Xt);
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


void material_damage::get_stress(const double *eps, double *stress,
				 const double *vars_old) const
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
	damage_law(eps, e_old, D_old, nullptr, &D, stress);

	for (int i = 0; i < 6; ++i)
		stress[i] *= (1 - D);
}


void material_damage::get_ctan(const double *eps, double *ctan,
			       const double *vars_old) const
{
	apply_perturbation(eps, ctan, vars_old);
}


bool material_damage::evolute(const double *eps, const double *vars_old,
			      double *vars_new) const
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

	bool non_linear = damage_law(eps, e_old, D_old, e_new, D_new, stress);

	return non_linear;
}


void material_damage::print() const
{
	cout << "Type : Damage" << endl;
	cout << "E = " << E << " nu = " << nu << " Xt = " << Xt << endl;
}

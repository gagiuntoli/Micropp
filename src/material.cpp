#include "material.hpp"

material_t * material_t::make_material(double *params, int type)
{
	/*
	 * This fabric creates the corresponding subclass according to the <type>
	 * and the <params> arguments
	 */
	switch (type) {
		case 0:
			return new material_elastic(params[0], params[1]);
			break;
		case 1:
			return new material_plastic(params[0], params[1], params[2], params[3]);
			break;
		case 2:
			return new material_damage(params[0], params[1], params[2]);
			break;
		default:
			break;
	}
}


material_t *material_t::make_material(material_t material)
{
	/*
	 * This fabric creates the corresponding subclass acorrding to the <type>
	 * data member of material.
	 */

	switch (material.type) {
		case 0:
			return new material_elastic(material.E, material.nu);
			break;
		case 1:
			return new material_plastic(material.E, material.nu, material.Ka, material.Sy);
			break;
		case 2:
			return new material_damage(material.E, material.nu, material.Xt);
			break;
		default:
			break;
	}
}


void get_dev_tensor(const double tensor[6], double tensor_dev[6])
{
	memcpy(tensor_dev, tensor, 6 * sizeof(double));
	for (int i = 0; i < 3; i++)
		tensor_dev[i] -= (1 / 3.0) * (tensor[0] + tensor[1] + tensor[2]);
}


// ELASTIC MATERIAL


void material_elastic::get_stress(const double *eps, double *stress, const double *history_params) const
{
	/* Elastic Material Law*/
	for (int i = 0; i < 3; ++i)
		stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = mu * eps[i];
}


void material_elastic::get_ctan(const double *eps, double *ctan, const double *history_params) const
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


bool material_elastic::evolute(const double *eps, const double *vars_old, double *vars_new) const
{
}


void material_elastic::print_n() const
{
	cout << "Type : Elastic" << endl;
	cout << "E = " << E << " nu = " << nu << endl;
}


// PLASTIC MATERIAL


bool material_plastic::plastic_law(const double eps[6],
				   const double *_eps_p_old,
				   const double *_alpha_old,
				   double *_dl,
				   double _normal[6],
				   double _s_trial[6]) const
{
	/* Calculates _dl, _normal and _s_trial
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


void material_plastic::get_ctan(const double *eps, double *ctan, const double *vars_old) const
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
			ctan[j * 6 + i] = (stress_1[j] - stress_0[j]) / D_EPS_CTAN;
	}
}


bool material_plastic::evolute(const double *eps, const double *vars_old, double *vars_new) const
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


void material_plastic::print_n() const
{
	cout << "Type : Plastic" << endl;
	cout << "E = " << E << " nu = " << nu << " Ka = " << Ka << " Sy = " << Sy << endl;
}


// DAMAGE MATERIAL


void material_damage::get_stress(const double *eps, double *stress, const double *history_params) const
{
	/* Elastic Material Law*/
	for (int i = 0; i < 3; ++i)
		stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = mu * eps[i];
}


void material_damage::get_ctan(const double *eps, double *ctan, const double *history_params) const
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


bool material_damage::evolute(const double *eps, const double *vars_old, double *vars_new) const
{
}


void material_damage::print_n() const
{
	cout << "Type : Plastic" << endl;
	cout << "E = " << E << " nu = " << nu << " Ka = " << Ka << " Sy = " << Sy << endl;
}

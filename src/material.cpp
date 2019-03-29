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


void material_elastic::get_stress(const double eps[6], double stress[6], const double *history_params) const
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


void material_elastic::print_n() const
{
	cout << "Type : Elastic" << endl;
	cout << "E = " << E << " nu = " << nu << endl;
}


void material_plastic::get_stress(const double eps[6], double stress[6], const double *history_params) const
{
	/* Elastic Material Law*/
	for (int i = 0; i < 3; ++i)
		stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * mu * eps[i];

	for (int i = 3; i < 6; ++i)
		stress[i] = mu * eps[i];
}


void material_plastic::get_ctan(const double *eps, double *ctan, const double *history_params) const
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


void material_plastic::print_n() const {
	cout << "Type : Plastic" << endl;
	cout << "E = " << E << " nu = " << nu << " Ka = " << Ka << " Sy = " << Sy << endl;
}


void material_damage::get_stress(const double eps[6], double stress[6], const double *history_params) const
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

void material_damage::print_n() const {
	cout << "Type : Plastic" << endl;
	cout << "E = " << E << " nu = " << nu << " Ka = " << Ka << " Sy = " << Sy << endl;
}

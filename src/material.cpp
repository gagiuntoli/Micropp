#include "material.hpp"

material_t * material_t::make_material(double * params, int type)
{
	switch (type) {
		case 0:
			return new material_elastic(params[0], params[1]);
			break;
		case 2:
			return new material_damage(params[0], params[1], params[2]);
			break;
		default:
			break;
	}
}

material_t *material_t::make_material(material_t *material)
{
	double params[4];
	switch (material->type) {
		case 0:
			params[0] = material->E;
			params[1] = material->nu;
			return new material_elastic(params[0], params[1]);
			break;
		case 1:
			params[0] = material->E;
			params[1] = material->nu;
			params[2] = material->Ka;
			params[3] = material->Sy;
			return new material_damage(params[0], params[1], params[2]);
			break;
		case 2:
			params[0] = material->E;
			params[1] = material->nu;
			params[2] = material->Xt;
			return new material_damage(params[0], params[1], params[2]);
			break;
		default:
			break;
	}
}


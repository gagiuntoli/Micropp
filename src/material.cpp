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


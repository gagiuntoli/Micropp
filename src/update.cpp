
#include "micropp.hpp"
#include "common.hpp"


/*
 * Evolutes the internal variables for the non-linear material models
 * Calculates the <f_trial_max> max value.
 */

template<int tdim>
bool micropp<tdim>::calc_vars_new(const double *u, const double *_vars_old,
				  double *_vars_new) const
{
	bool non_linear = false;

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex){

				const int e = glo_elem(ex, ey, ez);
				const material_t *material = get_material(e);

				for (int gp = 0; gp < npe; ++gp) {

					const double *vars_old = (_vars_old) ? &_vars_old[intvar_ix(e, gp, 0)] : nullptr;
					double *vars_new = &_vars_new[intvar_ix(e, gp, 0)];

					double eps[nvoi];
					get_strain(u, gp, eps, bmat, nx, ny, ex, ey, ez);

					non_linear |= material->evolute(eps, vars_old, vars_new);
				}
			}
		}
	}

	return non_linear;
}


template class micropp<3>;

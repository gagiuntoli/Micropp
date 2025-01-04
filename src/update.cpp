/*
 *  This source code is part of Micropp: a Finite Element library
 *  to solve composite materials micro-scale problems.
 *
 *  Copyright (C) - 2018
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

#include "common.hpp"
#include "micropp.hpp"

/*
 * Evolutes the internal variables for the non-linear material models
 * Calculates the <f_trial_max> max value.
 */

template <int tdim>
bool micropp<tdim>::calc_vars_new(const double *u, const double *_vars_old, double *_vars_new) const {
  bool non_linear = false;

  for (int ez = 0; ez < nez; ++ez) {
    for (int ey = 0; ey < ney; ++ey) {
      for (int ex = 0; ex < nex; ++ex) {
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

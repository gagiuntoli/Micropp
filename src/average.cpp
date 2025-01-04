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

template <int tdim>
void micropp<tdim>::calc_ave_strain(const double *u, double strain_ave[nvoi]) const {
  memset(strain_ave, 0, nvoi * sizeof(double));

  for (int ez = 0; ez < nez; ++ez) {  // 2D -> nez = 1
    for (int ey = 0; ey < ney; ++ey) {
      for (int ex = 0; ex < nex; ++ex) {
        double strain_aux[nvoi] = {0.0};

        for (int gp = 0; gp < npe; ++gp) {
          double strain_gp[nvoi];

          get_strain(u, gp, strain_gp, bmat, nx, ny, ex, ey, ez);
          for (int v = 0; v < nvoi; ++v) {
            strain_aux[v] += strain_gp[v] * wg;
          }
        }

        for (int v = 0; v < nvoi; v++) {
          strain_ave[v] += strain_aux[v];
        }
      }
    }
  }

  for (int v = 0; v < nvoi; v++) {
    strain_ave[v] /= vol_tot;
  }
}

template <int tdim>
void micropp<tdim>::calc_ave_stress(const double *u, double stress_ave[nvoi], const double *vars_old) const {
  memset(stress_ave, 0, nvoi * sizeof(double));

  for (int ez = 0; ez < nez; ++ez) {  // 2D -> nez = 1
    for (int ey = 0; ey < ney; ++ey) {
      for (int ex = 0; ex < nex; ++ex) {
        double stress_aux[nvoi] = {0.0};

        for (int gp = 0; gp < npe; ++gp) {
          double stress_gp[nvoi], strain_gp[nvoi];
          get_strain(u, gp, strain_gp, bmat, nx, ny, ex, ey, ez);
          get_stress(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);
          for (int v = 0; v < nvoi; ++v) {
            stress_aux[v] += stress_gp[v] * wg;
          }
        }
        for (int v = 0; v < nvoi; ++v) {
          stress_ave[v] += stress_aux[v];
        }
      }
    }
  }

  for (int v = 0; v < nvoi; ++v) stress_ave[v] /= vol_tot;
}

template <int tdim>
void micropp<tdim>::calc_fields(double *u, double *vars_old) {
  for (int ez = 0; ez < nez; ++ez) {  // 2D -> nez = 1
    for (int ey = 0; ey < ney; ++ey) {
      for (int ex = 0; ex < nex; ++ex) {
        double eps_a[nvoi] = {0.0};
        double sig_a[nvoi] = {0.0};

        for (int gp = 0; gp < npe; ++gp) {
          double stress_gp[nvoi], strain_gp[nvoi];

          get_strain(u, gp, strain_gp, bmat, nx, ny, ex, ey, ez);
          get_stress(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);

          for (int v = 0; v < nvoi; ++v) {
            eps_a[v] += strain_gp[v] * wg;
            sig_a[v] += stress_gp[v] * wg;
          }
        }

        const int e = glo_elem(ex, ey, ez);
        for (int v = 0; v < nvoi; ++v) {
          elem_strain[e * nvoi + v] = eps_a[v] * ivol;
          elem_stress[e * nvoi + v] = sig_a[v] * ivol;
        }
      }
    }
  }
}

template class micropp<3>;

/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *						             Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <cmath>

#include "instrument.hpp"
#include "micropp.hpp"

template <int tdim>
void micropp<tdim>::set_strain(const int gp_id, const double *strain) {
  assert(gp_id >= 0);
  assert(gp_id < ngp);
  memcpy(gp_list[gp_id].strain, strain, nvoi * sizeof(double));
}

template <int tdim>
void micropp<tdim>::get_stress(const int gp_id, double *stress) const {
  assert(gp_id >= 0);
  assert(gp_id < ngp);
  memcpy(stress, gp_list[gp_id].stress, nvoi * sizeof(double));
}

template <int tdim>
void micropp<tdim>::get_ctan(const int gp_id, double *ctan) const {
  assert(gp_id >= 0);
  assert(gp_id < ngp);
  memcpy(ctan, gp_list[gp_id].ctan, nvoi * nvoi * sizeof(double));
}

template <int tdim>
void micropp<tdim>::homogenize_linear() {
  INST_START;

#pragma omp parallel for schedule(dynamic, 1)
  for (int igp = 0; igp < ngp; ++igp) {
    gp_t<tdim> *gp_ptr = &gp_list[igp];

    /*
     * Computational cheap calculation
     * stress = ctan_lin * strain
     */

    homogenize_linear(gp_ptr);
  }
}

template <int tdim>
void micropp<tdim>::homogenize() {
  INST_START;

#pragma omp parallel for schedule(dynamic, 1)
  for (int igp = 0; igp < ngp; ++igp) {
    gp_t<tdim> *gp_ptr = &gp_list[igp];

    if (gp_ptr->coupling == FE_LINEAR || gp_ptr->coupling == MIX_RULE_CHAMIS) {
      /*
       * Computational cheap calculation
       * stress = ctan_lin * strain
       *
       * All mixture rules are linear in Micropp
       * so the homogenization of the stress tensor
       * is this simple and cheap procedure.
       */

      homogenize_linear(gp_ptr);

    } else if (gp_ptr->coupling == FE_ONE_WAY) {
      homogenize_fe_one_way(gp_ptr);

    } else if (gp_ptr->coupling == FE_FULL) {
      homogenize_fe_full(gp_ptr);
    }
  }

  if (write_log_flag) {
    write_log();
  }
}

template <int tdim>
void micropp<tdim>::homogenize_linear(gp_t<tdim> *gp_ptr) {
  memset(gp_ptr->stress, 0.0, nvoi * sizeof(double));
  for (int i = 0; i < nvoi; ++i) {
    for (int j = 0; j < nvoi; ++j) {
      gp_ptr->stress[i] += gp_ptr->ctan[i * nvoi + j] * gp_ptr->strain[j];
    }
  }
}

template <int tdim>
void micropp<tdim>::homogenize_fe_one_way(gp_t<tdim> *gp_ptr) {
  ell_matrix A;  // Jacobian
  const int ns[3] = {nx, ny, nz};
  ell_init(&A, dim, dim, ns, CG_ABS_TOL, CG_REL_TOL, CG_MAX_ITS);
  double *b = (double *)calloc(nndim, sizeof(double));
  double *du = (double *)calloc(nndim, sizeof(double));
  double *u = (double *)calloc(nndim, sizeof(double));
  double *vars_new_aux = (double *)calloc(nvars, sizeof(double));

  double *vars_new = (gp_ptr->allocated) ? gp_ptr->vars_k : vars_new_aux;

  gp_ptr->cost = 0;
  gp_ptr->subiterated = false;

  // SIGMA 1 Newton-Raphson
  memcpy(u, gp_ptr->u_n, nndim * sizeof(double));

  newton_t newton = newton_raphson(&A, b, u, du, gp_ptr->strain, gp_ptr->vars_n);

  memcpy(gp_ptr->u_k, u, nndim * sizeof(double));
  gp_ptr->cost += newton.solver_its;
  gp_ptr->converged = newton.converged;

  /*
   * In case it has not converged do the sub-iterations
   */
  if (gp_ptr->converged == false && subiterations == true) {
    double eps_sub[nvoi], deps_sub[nvoi];
    memcpy(u, gp_ptr->u_n, nndim * sizeof(double));
    memcpy(eps_sub, gp_ptr->strain_old, nvoi * sizeof(double));
    gp_ptr->subiterated = true;

    for (int i = 0; i < nvoi; ++i) deps_sub[i] = (gp_ptr->strain[i] - gp_ptr->strain_old[i]) / nsubiterations;

    for (int its = 0; its < nsubiterations; ++its) {
      for (int j = 0; j < nvoi; ++j) {
        eps_sub[j] += deps_sub[j];
      }

      newton = newton_raphson(&A, b, u, du, eps_sub, gp_ptr->vars_n);
      gp_ptr->cost += newton.solver_its;
    }

    gp_ptr->converged = newton.converged;
    memcpy(gp_ptr->u_k, u, nndim * sizeof(double));
  }

  if (lin_stress) {
    memset(gp_ptr->stress, 0.0, nvoi * sizeof(double));
    for (int i = 0; i < nvoi; ++i) {
      for (int j = 0; j < nvoi; ++j) {
        gp_ptr->stress[i] += gp_ptr->ctan[i * nvoi + j] * gp_ptr->strain[j];
      }
    }

  } else {
    calc_ave_stress(gp_ptr->u_k, gp_ptr->stress, gp_ptr->vars_n);
  }

  // Updates <vars_new>
  bool non_linear = calc_vars_new(gp_ptr->u_k, gp_ptr->vars_n, vars_new);

  if (non_linear == true) {
    if (gp_ptr->allocated == false) {
      gp_ptr->allocate();
      memcpy(gp_ptr->vars_k, vars_new, nvars * sizeof(double));
    }
  }

  ell_free(&A);
  free(b);
  free(u);
  free(du);
  free(vars_new_aux);
}

template <int tdim>
void micropp<tdim>::homogenize_fe_full(gp_t<tdim> *gp_ptr) {
  ell_matrix A;  // Jacobian
  const int ns[3] = {nx, ny, nz};
  ell_init(&A, dim, dim, ns, CG_ABS_TOL, CG_REL_TOL, CG_MAX_ITS);
  double *b = (double *)calloc(nndim, sizeof(double));
  double *du = (double *)calloc(nndim, sizeof(double));
  double *u = (double *)calloc(nndim, sizeof(double));
  double *vars_new_aux = (double *)calloc(nvars, sizeof(double));

  double *vars_new = (gp_ptr->allocated) ? gp_ptr->vars_k : vars_new_aux;

  gp_ptr->cost = 0;
  gp_ptr->subiterated = false;

  // SIGMA 1 Newton-Raphson
  memcpy(u, gp_ptr->u_n, nndim * sizeof(double));

  newton_t newton = newton_raphson(&A, b, u, du, gp_ptr->strain, gp_ptr->vars_n);

  memcpy(gp_ptr->u_k, u, nndim * sizeof(double));
  gp_ptr->cost += newton.solver_its;
  gp_ptr->converged = newton.converged;

  /*
   * In case it has not converged do the sub-iterations
   */
  if (gp_ptr->converged == false && subiterations == true) {
    double eps_sub[nvoi], deps_sub[nvoi];
    memcpy(u, gp_ptr->u_n, nndim * sizeof(double));
    memcpy(eps_sub, gp_ptr->strain_old, nvoi * sizeof(double));
    gp_ptr->subiterated = true;

    for (int i = 0; i < nvoi; ++i) deps_sub[i] = (gp_ptr->strain[i] - gp_ptr->strain_old[i]) / nsubiterations;

    for (int its = 0; its < nsubiterations; ++its) {
      for (int j = 0; j < nvoi; ++j) eps_sub[j] += deps_sub[j];

      newton = newton_raphson(&A, b, u, du, eps_sub, gp_ptr->vars_n);
      gp_ptr->cost += newton.solver_its;
    }

    gp_ptr->converged = newton.converged;
    memcpy(gp_ptr->u_k, u, nndim * sizeof(double));
  }

  if (lin_stress) {
    memset(gp_ptr->stress, 0.0, nvoi * sizeof(double));
    for (int i = 0; i < nvoi; ++i) {
      for (int j = 0; j < nvoi; ++j) {
        gp_ptr->stress[i] += gp_ptr->ctan[i * nvoi + j] * gp_ptr->strain[j];
      }
    }

  } else {
    calc_ave_stress(gp_ptr->u_k, gp_ptr->stress, gp_ptr->vars_n);
  }

  // Updates <vars_new>
  bool non_linear = calc_vars_new(gp_ptr->u_k, gp_ptr->vars_n, vars_new);

  if (non_linear == true) {
    if (gp_ptr->allocated == false) {
      gp_ptr->allocate();
      memcpy(gp_ptr->vars_k, vars_new, nvars * sizeof(double));
    }
  }

  if (gp_ptr->allocated) {
    // CTAN 3/6 Newton-Raphsons in 2D/3D
    double eps_1[6], sig_0[6], sig_1[6];

    memcpy(u, gp_ptr->u_k, nndim * sizeof(double));
    memcpy(sig_0, gp_ptr->stress, nvoi * sizeof(double));

    for (int i = 0; i < nvoi; ++i) {
      memcpy(eps_1, gp_ptr->strain, nvoi * sizeof(double));
      eps_1[i] += D_EPS_CTAN_AVE;

      newton_raphson(&A, b, u, du, eps_1, gp_ptr->vars_n);

      gp_ptr->cost += newton.solver_its;

      calc_ave_stress(u, sig_1, gp_ptr->vars_n);

      for (int v = 0; v < nvoi; ++v) gp_ptr->ctan[v * nvoi + i] = (sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;
    }
  }

  ell_free(&A);
  free(b);
  free(u);
  free(du);
  free(vars_new_aux);
}

template <int tdim>
void micropp<tdim>::update_vars() {
  for (int igp = 0; igp < ngp; ++igp) gp_list[igp].update_vars();
}

template class micropp<3>;

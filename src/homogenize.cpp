/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
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

#include <cmath>
#include <cassert>

#include "instrument.hpp"
#include "micro.hpp"


template <int tdim>
void micropp<tdim>::set_macro_strain(const int gp_id,
                                     const double *macro_strain)
{
    INST_START;

    assert(gp_id >= 0);
    assert(gp_id < ngp);
    memcpy(gp_list[gp_id].macro_strain, macro_strain, nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_stress(const int gp_id,
                                     double *macro_stress) const
{
    INST_START;

    assert(gp_id >= 0);
    assert(gp_id < ngp);
    memcpy(macro_stress, gp_list[gp_id].macro_stress, nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_ctan(const int gp_id, double *macro_ctan) const
{
    INST_START;

    assert(gp_id >= 0);
    assert(gp_id < ngp);
    memcpy(macro_ctan, gp_list[gp_id].macro_ctan, nvoi * nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::homogenize()
{
    INST_START;

    int newton_its, solver_its[NR_MAX_ITS];
    double newton_err[NR_MAX_ITS], solver_err[NR_MAX_ITS];
    bool nl_flag;

    for (int igp = 0; igp < ngp; ++igp) {
        gp_t<tdim> * const gp_ptr = &gp_list[igp];

        gp_ptr->sigma_cost = 0;

        if (!gp_ptr->allocated) {
            vars_old = vars_old_aux;
            vars_new = vars_new_aux;
            memset(vars_old, 0, num_int_vars * sizeof(double));
        } else {
            vars_old = gp_ptr->int_vars_n;
            vars_new = gp_ptr->int_vars_k;
        }

        // SIGMA 1 Newton-Raphson
        memcpy(gp_ptr->u_k, gp_ptr->u_n, nndim * sizeof(double));

        newton_its = newton_raphson(gp_ptr->macro_strain, gp_ptr->u_k,
                                    newton_err, solver_its, solver_err);

        gp_ptr->sigma_newton_its = newton_its;
        memcpy(gp_ptr->sigma_newton_err, newton_err, NR_MAX_ITS * sizeof(double));
        memcpy(gp_ptr->sigma_solver_its, solver_its,
               NR_MAX_ITS * sizeof(int));
        memcpy(gp_ptr->sigma_solver_err, solver_err,
               NR_MAX_ITS * sizeof(double));

        for (int i = 0; i < newton_its; ++i)
            gp_ptr->sigma_cost += solver_its[i];

        calc_ave_stress(gp_ptr->u_k, gp_ptr->macro_stress);

        nl_flag = calc_vars_new(gp_ptr->u_k);

        if (nl_flag) {
            if (!gp_ptr->allocated) {
                gp_ptr->allocate(num_int_vars);
                memcpy(gp_ptr->int_vars_k, vars_new,
                       num_int_vars * sizeof(double));
            }
        }

        if (gp_ptr->allocated) {

            // CTAN 3/6 Newton-Raphsons in 2D/3D
            double eps_1[6], sig_0[6], sig_1[6];

            memcpy(u_aux, gp_ptr->u_k, nndim * sizeof(double));
            memcpy(sig_0, gp_ptr->macro_stress, nvoi * sizeof(double));

            for (int i = 0; i < nvoi; ++i) {

                memcpy(eps_1, gp_ptr->macro_strain, nvoi * sizeof(double));
                eps_1[i] += D_EPS_CTAN_AVE;

                newton_its = newton_raphson(eps_1, u_aux, newton_err,
                                            solver_its, solver_err);

                gp_ptr->nr_its[i + 1] = newton_its;
                gp_ptr->nr_err[i + 1] = newton_err[0];

                calc_ave_stress(u_aux, sig_1);

                for (int v = 0; v < nvoi; ++v)
                    gp_ptr->macro_ctan[v * nvoi + i] =
                        (sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;

            }
        }
    }
}


template <int tdim>
void micropp<tdim>::update_vars()
{
    INST_START;

    for (int igp = 0; igp < ngp; ++igp)
        gp_list[igp].update_vars();
}


template class micropp<2>;
template class micropp<3>;

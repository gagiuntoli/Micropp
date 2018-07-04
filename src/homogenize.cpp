/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
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

#include "micro.h"

void micropp_t::calc_ctan_lin()
{
	double sig_1[6], eps_1[6];
	double d_eps = 1.0e-8;

	for (int i=0; i<nvoi; i++)
	{
		for (int v=0; v<nvoi; v++)
  			eps_1[v] = 0.0;

		eps_1[i] += d_eps;

		int nr_its;
		bool nl_flag;
        double nr_err;
		set_displ(eps_1);
		newton_raphson(&nl_flag, &nr_its, &nr_err);
		calc_ave_stress(sig_1);

		for (int v = 0; v < nvoi; ++v)
  			ctan_lin[v*nvoi + i] = sig_1[v] / d_eps;
	}
}

bool micropp_t::is_linear (const double *macro_strain)
{
	double macro_stress[6];
	for (int i = 0; i < nvoi; ++i) {
		macro_stress[i] = 0.0;
		for (int j = 0; j < nvoi; ++j)
  			macro_stress[i] += ctan_lin[i*nvoi + j] * macro_strain[j];
	}

	double inv = get_inv_1(macro_stress);
	if (fabs(inv) > inv_max)
 		inv_max = fabs(inv);

	return (fabs(inv) < inv_tol) ? true : false;
}

double micropp_t::get_inv_1(const double *tensor)
{
	if (dim == 2)
		return tensor[0] + tensor[1];
	if (dim == 3)
		return tensor[0] + tensor[1] + tensor[2];
}

double micropp_t::get_inv_2(const double *tensor)
{
	if (dim == 3)
		return \
			tensor[0]*tensor[1] + tensor[0]*tensor[2] + tensor[1]*tensor[2] + \
			tensor[3]*tensor[3] + tensor[4]*tensor[4] + tensor[5]*tensor[5];
}

void micropp_t::set_macro_strain(const int gp_id, const double *macro_strain)
{
	list<gp_t>::iterator gp;
	for(gp = gauss_list.begin(); gp != gauss_list.end(); ++gp) {
		if (gp->id == gp_id) {
  			for (int i = 0; i < nvoi; ++i)
				gp->macro_strain[i] = macro_strain[i];
  			break;
		}
	}
	if (gp == gauss_list.end()) {
		gp_t gp_n;
		gp_n.id = gp_id;
		gp_n.int_vars_n = NULL;
		gp_n.int_vars_k = NULL;
		for (int i=0; i<nvoi; i++)
  			gp_n.macro_strain[i] = macro_strain[i];
		gp_n.inv_max = -1.0e10;
		gauss_list.push_back(gp_n);
	}
}

void micropp_t::get_macro_stress(const int gp_id, double *macro_stress)
{
	for (auto const& gp : gauss_list)
		if (gp.id == gp_id) {
  			for (int i=0; i<nvoi; i++)
				macro_stress[i] = gp.macro_stress[i];
  			break;
		}
}

void micropp_t::get_macro_ctan(const int gp_id, double *macro_ctan)
{
	for (auto const& gp : gauss_list)
		if (gp.id == gp_id) {
  			for (int i = 0; i < (nvoi*nvoi); ++i)
				macro_ctan[i] = gp.macro_ctan[i];
  			break;
		}
}

void micropp_t::homogenize()
{
	for (auto& gp : gauss_list) {

		if (gp.int_vars_n == NULL)
  			for (int i = 0; i < num_int_vars; ++i)
				vars_old[i] = 0.0;
		else
  			for (int i = 0; i < num_int_vars; ++i)
				vars_old[i] = gp.int_vars_n[i];

		inv_max = -1.0e10;

		if ((is_linear(gp.macro_strain) == true) && (gp.int_vars_n == NULL)) {

  			// S = CL : E
  			for (int i = 0; i < nvoi; ++i) {
				gp.macro_stress[i] = 0.0;
				for (int j = 0; j < nvoi; ++j)
  					gp.macro_stress[i] += ctan_lin[i*nvoi + j] * gp.macro_strain[j];
  			}
  			// C = CL
  			for (int i = 0; i < nvoi; ++i)
				for (int j = 0; j < nvoi; ++j)
  					gp.macro_ctan[i*nvoi + j] = ctan_lin[i*nvoi + j];

  			for (int i = 0; i < (1+nvoi); ++i) {
				gp.nr_its[i] = 0;
				gp.nr_err[i] = 0.0;
  			}

		} else {

		 	// SIGMA
			int nr_its;
  			bool nl_flag;
        	double nr_err;
  			set_displ(gp.macro_strain);
  			newton_raphson(&nl_flag, &nr_its, &nr_err); 
  			calc_ave_stress(gp.macro_stress);

			if (nl_flag == true)
			{
  				if(gp.int_vars_n == NULL) {
					gp.int_vars_k = (double*)malloc(num_int_vars*sizeof(double));
					gp.int_vars_n = (double*)malloc(num_int_vars*sizeof(double));
  					for (int i = 0; i < num_int_vars; ++i)
						gp.int_vars_n[i] = 0.0;
				}
  				for (int i = 0; i < num_int_vars; ++i)
					gp.int_vars_k[i] = vars_new[i];
			}

  			gp.nr_its[0] = nr_its;
  			gp.nr_err[0] = nr_err;

  			// CTAN
  			double eps_1[6], sig_0[6], sig_1[6], dEps = 1.0e-8;
			for (int v = 0; v < nvoi; ++v)
  				sig_0[v] = gp.macro_stress[v];

  			for (int i = 0; i < nvoi; ++i)
 			{
				for (int v = 0; v < nvoi; ++v)
  					eps_1[v] = gp.macro_strain[v];
				eps_1[i] += dEps;

				set_displ(eps_1);
				newton_raphson(&nl_flag, &nr_its, &nr_err);
				calc_ave_stress(sig_1);
				for (int v = 0; v < nvoi; ++v)
  					gp.macro_ctan[v*nvoi + i] = (sig_1[v] - sig_0[v]) / dEps;

				gp.nr_its[1+i] = nr_its;
				gp.nr_err[1+i] = nr_err;
  			}
		}
		gp.inv_max = inv_max;
	}
}

void micropp_t::update_vars()
{
	for (auto const& gp : gauss_list)
		if(gp.int_vars_n != NULL)
  			for (int i = 0; i < num_int_vars; ++i)
				gp.int_vars_n[i] = gp.int_vars_k[i];
}

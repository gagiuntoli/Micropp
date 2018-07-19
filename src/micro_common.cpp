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

#include "micro.hpp"

template<int tdim>
void micropp<tdim>::initialize(const double *_micro_params,
                               const int *_mat_types,
                               const double *_params)
{
	INST_CONSTRUCT; // Initialize the Intrumentation

	b = (double *) malloc(nn * dim * sizeof(double));
	du = (double *) malloc(nn * dim * sizeof(double));
	u = (double *) malloc(nn * dim * sizeof(double));

	gp_list = new gp_t[ngp]();

	elem_type = (int *) malloc(nelem * sizeof(int));
	elem_stress = (double *) malloc(nelem * nvoi * sizeof(double));
	elem_strain = (double *) malloc(nelem * nvoi * sizeof(double));
	vars_old_aux = (double *) calloc(num_int_vars, sizeof(double));
	vars_new_aux = (double *) malloc(num_int_vars* sizeof(double));
	vars_old = vars_old_aux;
	vars_new = vars_new_aux;

	output_files_header = false;
	assert( b && du && u && elem_stress && elem_strain &&
	        elem_type && vars_old && vars_new );

	int nParams;
	if (micro_type == 0) {
		// mat 1 = matrix
		// mat 2 = sphere
		numMaterials = 2;
		nParams = 5;
	} else if (micro_type == 1) {
		// mat 1 = layer 1
		// mat 2 = layer 2
		numMaterials = 2;
		nParams = 5;
	}

	for (int i = 0; i < nParams; i++)
		micro_params[i] = _micro_params[i];

	ofstream file;
	file.open("micropp_materials.dat");
	file << scientific;
	for (int i = 0; i < numMaterials; i++) {
		material_list[i].E  = _params[i * MAX_MAT_PARAM + 0];
		material_list[i].nu = _params[i * MAX_MAT_PARAM + 1];
		material_list[i].Sy = _params[i * MAX_MAT_PARAM + 2];
		material_list[i].Ka = _params[i * MAX_MAT_PARAM + 3];

		material_list[i].k = material_list[i].E / (3 * (1 - 2 * material_list[i].nu));
		material_list[i].mu = material_list[i].E / (2 * (1 + material_list[i].nu));
		material_list[i].lambda = (material_list[i].nu * material_list[i].E) /
			((1 + material_list[i].nu) * (1 - 2 * material_list[i].nu));	// lambda

		if (_mat_types[i] == 0) {
			// lineal
			material_list[i].plasticity = false;
			material_list[i].damage = false;
		} else if (_mat_types[i] == 1) {
			// con plasticidad
			material_list[i].plasticity = true;
			material_list[i].damage = false;
		} else if (_mat_types[i] == 2) {
			// con daño
			material_list[i].plasticity = false;
			material_list[i].damage = true;
		}
		file << setw(14) << material_list[i].E << " " << material_list[i].nu
		     << " " << material_list[i].Sy << " " << material_list[i].Ka << endl;
	}
	file.close();

	solver.max_its = CG_MAX_ITS;
	solver.min_tol = CG_MAX_TOL;
	solver.k = (double *) malloc(nn * dim * sizeof(double));
	solver.r = (double *) malloc(nn * dim * sizeof(double));
	solver.z = (double *) malloc(nn * dim * sizeof(double));
	solver.p = (double *) malloc(nn * dim * sizeof(double));
	solver.q = (double *) malloc(nn * dim * sizeof(double));

	file.open("micropp_convergence.dat");
	file.close();
	file.open("micropp_eps_sig_ctan.dat");
	file.close();
	file.open("micropp_int_vars_n.dat");
	file.close();
}

template <int tdim>
micropp<tdim>::~micropp()
{
	INST_DESTRUCT;

	ell_free(&A);

	free(b);
	free(du);
	free(u);
	free(elem_stress);
	free(elem_strain);
	free(elem_type);
	free(vars_old_aux);
	free(vars_new_aux);

	delete [] gp_list;

	free(solver.k);
	free(solver.r);
	free(solver.z);
	free(solver.p);
	free(solver.q);
}


template <int tdim>
int micropp<tdim>::get_nl_flag(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].allocated;
}


template <int tdim>
void micropp<tdim>::calc_ctan_lin()
{
	double sig_1[6], eps_1[6];
	double d_eps = 1.0e-8;

	for (int i = 0; i < nvoi; i++) {
		for (int v = 0; v < nvoi; v++)
			eps_1[v] = 0.0;

		eps_1[i] += d_eps;

		int nr_its;
		bool nl_flag;
		double nr_err;
		set_displ(eps_1);
		newton_raphson(&nl_flag, &nr_its, &nr_err);
		calc_ave_stress(sig_1);

		for (int v = 0; v < nvoi; ++v)
			ctan_lin[v * nvoi + i] = sig_1[v] / d_eps;
	}
}


template <int tdim>
bool micropp<tdim>::is_linear(const double *macro_strain)
{
	double macro_stress[6];
	for (int i = 0; i < nvoi; ++i) {
		macro_stress[i] = 0.0;
		for (int j = 0; j < nvoi; ++j)
			macro_stress[i] += ctan_lin[i * nvoi + j] * macro_strain[j];
	}

	const double inv = get_inv_1(macro_stress);
	if (fabs(inv) > inv_max)
		inv_max = fabs(inv);

	return (fabs(inv) < inv_tol);
}


template <int tdim>
double micropp<tdim>::get_inv_1(const double *tensor) const
{
	assert(dim == 2 || dim == 3);

	const double ret = tensor[0] + tensor[1];

	if (dim == 2)
		return ret;

	return ret + tensor[2];
}


template <int tdim>
material_t micropp<tdim>::get_material(const int e) const
{
	int mat_num;
	if (micro_type == 0)
		if (elem_type[e] == 0)
			mat_num = 0;
		else
			mat_num = 1;

	else if (micro_type == 1)
		if (elem_type[e] == 0)
			mat_num = 0;
		else
			mat_num = 1;

	return material_list[mat_num];
}


template <int tdim>
void micropp<tdim>::print_info() const
{
	printf("micropp%d\n", dim);
	printf("ngp %d n = [%d, %d, %d] => nn = %d\n", ngp, nx, ny, nz, nn);
	printf("l = [%lf, %lf, %lf]; width = %lf\n", lx, ly, lz, width);
}

// Explicit instantiation
template class micropp<2>;
template class micropp<3>;

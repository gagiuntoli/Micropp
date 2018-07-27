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
                               const material_t *_materials)
{
	INST_CONSTRUCT; // Initialize the Intrumentation

	gp_list = new gp_t[ngp]();
	for (int gp = 0; gp < ngp; gp++) {
		gp_list[gp].u_n = (double *) calloc(nn * dim, sizeof(double));
		gp_list[gp].u_k = (double *) malloc(nn * dim * sizeof(double));
	}

	b = (double *) malloc(nn * dim * sizeof(double));
	du = (double *) malloc(nn * dim * sizeof(double));
	u_aux = (double *) malloc(nn * dim * sizeof(double));

	elem_type = (int *) malloc(nelem * sizeof(int));
	elem_stress = (double *) malloc(nelem * nvoi * sizeof(double));
	elem_strain = (double *) malloc(nelem * nvoi * sizeof(double));
	vars_old_aux = (double *) calloc(num_int_vars, sizeof(double));
	vars_new_aux = (double *) malloc(num_int_vars* sizeof(double));

	assert(b && du && u_aux && elem_stress && elem_strain &&
			elem_type && vars_old_aux && vars_new_aux);

	output_files_header = false;

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
		material_list[i] = _materials[i];

		file << setw(14)
		     << material_list[i].E << " " << material_list[i].nu << " "
		     << material_list[i].Sy << " " << material_list[i].Ka << endl;
	}

	for (int ez = 0; ez < nez; ez++) {
		for (int ey = 0; ey < ney; ey++) {
			for (int ex = 0; ex < nex; ex++) {
				const int e_i = glo_elem(ex, ey, ez);
				elem_type[e_i] = get_elem_type(ex, ey, ez);
			}
		}
	}

	const int ns[3] = { nx, ny, nz };
	const int nfield = dim;
	ell_init(&A, nfield, dim, ns, CG_MIN_ERR, CG_MAX_ITS);

	calc_ctan_lin();

	file.close();

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
	free(u_aux);
	free(elem_stress);
	free(elem_strain);
	free(elem_type);
	free(vars_old_aux);
	free(vars_new_aux);

	delete [] gp_list;
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
	u = u_aux;
	vars_old = vars_old_aux;

	for (int i = 0; i < nvoi; ++i) {

        double eps_1[nvoi] = { 0.0 };
		eps_1[i] += D_EPS_CTAN_AVE;

		double nr_err;
		set_displ_bc(eps_1);
		newton_raphson(&nr_err);

		double sig_1[6];
		calc_ave_stress(sig_1);

		for (int v = 0; v < nvoi; ++v)
			ctan_lin[v * nvoi + i] = sig_1[v] / D_EPS_CTAN_AVE;
	}
}


template <int tdim>
bool micropp<tdim>::is_linear(const double *macro_strain)
{
	double macro_stress[6] = { 0.0 };
	for (int i = 0; i < nvoi; ++i)
		for (int j = 0; j < nvoi; ++j)
			macro_stress[i] += ctan_lin[i * nvoi + j] * macro_strain[j];

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
void micropp<tdim>::get_elem_nodes(int n[8], int ex, int ey, int ez) const
{
	const int nxny = ny * nx;
	const int n0 = ez * nxny + ey * nx + ex;
	n[0] = n0;
	n[1] = n0 + 1;
	n[2] = n0 + nx + 1;
	n[3] = n0 + nx;
	n[4] = n[0] + nxny;
	n[5] = n[1] + nxny;
	n[6] = n[2] + nxny;
	n[7] = n[3] + nxny;
}


template <int tdim>
void micropp<tdim>::get_elem_displ(const double *u,
		double *elem_disp, int ex, int ey, int ez) const
{
	int n[8] ;
	get_elem_nodes(n, ex, ey, ez);

	for (int i = 0 ; i < npe; ++i)
		for (int d = 0; d < dim; ++d)
			elem_disp[i * dim + d] = u[n[i] * dim + d];
}


template <int tdim>
void micropp<tdim>::get_strain(int gp, double *strain_gp,
		int ex, int ey, int ez) const
{
	double elem_disp[3 * 8];
	get_elem_displ(u, elem_disp, ex, ey, ez);

	double bmat[6][3 * 8];
	calc_bmat(gp, bmat);

	memset(strain_gp, 0, nvoi * sizeof(double));
	for (int v = 0; v < nvoi; ++v)
		for (int i = 0; i < npe * dim; i++)
			strain_gp[v] += bmat[v][i] * elem_disp[i];
}


template <int tdim>
void micropp<tdim>::print_info() const
{
	printf("micropp%d\n", dim);
	printf("ngp %d n = [%d, %d, %d] => nn = %d\n", ngp, nx, ny, nz, nn);
	printf("l = [%lf, %lf, %lf]; width = %lf\n", lx, ly, lz, width);
}


template <int tdim>
void micropp<tdim>::get_stress(int gp, const double eps[nvoi],
		double stress_gp[nvoi], int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t material = get_material(e);
	const double mu = material.mu;

	if (material.plasticity == true) {

		const double *eps_p_old = &vars_old[intvar_ix(e, gp, 0)];
		const double alpha_old = vars_old[intvar_ix(e, gp, 6)];

		plastic_get_stress(&material, eps, eps_p_old, alpha_old, stress_gp);

	} else {

		isolin_get_stress(&material, eps, stress_gp);
	}

}

// Explicit instantiation
template class micropp<2>;
template class micropp<3>;

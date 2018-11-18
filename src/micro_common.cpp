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
micropp<tdim>::micropp(const int _ngp, const int size[3], const int _micro_type,
		       const double _micro_params[4],
		       const material_t *_materials):
	ngp(_ngp),
	nx(size[0]), ny(size[1]),
	nz((tdim == 3) ? size[2] : 1),

	nn(nx * ny * nz),
	nndim(nn * dim),

	nex(nx - 1), ney(ny - 1),
	nez((tdim == 3) ? (nz - 1) : 1),

	nelem(nex * ney * nez),
	lx(_micro_params[0]), ly(_micro_params[1]),
	lz((tdim == 3) ? _micro_params[2] : 0.0),
	dx(lx / nex), dy(ly / ney), dz((tdim == 3) ? lz / nez : 0.0),

	special_param(_micro_params[3]),

	wg(((tdim == 3) ? dx * dy * dz : dx * dy) / npe),
	vol_tot((tdim == 3) ? lx * ly * lz : lx * ly),
	ivol(1.0 / (wg * npe)),
	micro_type(_micro_type), num_int_vars(nelem * npe * NUM_VAR_GP)
{
	INST_CONSTRUCT; // Initialize the Intrumentation

	gp_list = new gp_t<tdim>[ngp]();
	for (int gp = 0; gp < ngp; ++gp) {
		gp_list[gp].u_n = (double *) calloc(nndim, sizeof(double));
		gp_list[gp].u_k = (double *) calloc(nndim, sizeof(double));
	}

	b = (double *) calloc(nndim, sizeof(double));
	du = (double *) calloc(nndim, sizeof(double));
	u_aux = (double *) calloc(nndim, sizeof(double));

	elem_type = (int *) calloc(nelem, sizeof(int));
	elem_stress = (double *) calloc(nelem * nvoi, sizeof(double));
	elem_strain = (double *) calloc(nelem * nvoi, sizeof(double));
	vars_new_aux = (double *) calloc(num_int_vars, sizeof(double));

	assert(b && du && u_aux && elem_stress
	       && elem_strain && elem_type && vars_new_aux);

	int nParams = 4;
	numMaterials = 2;

	for (int i = 0; i < nParams; i++)
		micro_params[i] = _micro_params[i];

	for (int i = 0; i < numMaterials; ++i)
		material_list[i] = _materials[i];

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {
				const int e_i = glo_elem(ex, ey, ez);
				elem_type[e_i] = get_elem_type(ex, ey, ez);
			}
		}
	}

	const int ns[3] = { nx, ny, nz };
	const int nfield = dim;
	ell_init(&A, nfield, dim, ns, CG_MIN_ERR, CG_MAX_ITS);
	ell_init(&A0, nfield, dim, ns, CG_MIN_ERR, CG_MAX_ITS);
	assembly_mat(&A0, u_aux, NULL);

	calc_ctan_lin();

	for (int gp = 0; gp < ngp; ++gp)
		gp_list[gp].macro_ctan = ctan_lin;

}


template <int tdim>
micropp<tdim>::~micropp()
{
	INST_DESTRUCT;

	ell_free(&A);
	ell_free(&A0);

	free(b);
	free(du);
	free(u_aux);
	free(elem_stress);
	free(elem_strain);
	free(elem_type);
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
void micropp<tdim>::get_sigma_solver_its(int gp_id,
					 int sigma_solver_its[NR_MAX_ITS])
	const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	memcpy(sigma_solver_its, gp_list[gp_id].sigma_solver_its,
	       NR_MAX_ITS * sizeof(int));
}


	template <int tdim>
void micropp<tdim>::get_sigma_solver_err(int gp_id,
					 double sigma_solver_err[NR_MAX_ITS])
	const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	memcpy(sigma_solver_err, gp_list[gp_id].sigma_solver_err,
	       NR_MAX_ITS * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_sigma_newton_err(int gp_id,
					 double sigma_newton_err[NR_MAX_ITS])
	const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	memcpy(sigma_newton_err, gp_list[gp_id].sigma_newton_err,
	       NR_MAX_ITS * sizeof(double));
}


template <int tdim>
int micropp<tdim>::get_sigma_newton_its(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].sigma_newton_its;
}


template <int tdim>
int micropp<tdim>::get_sigma_cost(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].sigma_cost;
}


template <int tdim>
void micropp<tdim>::calc_ctan_lin()
{
	double sig_1[6];

	for (int i = 0; i < nvoi; ++i) {

		double eps_1[nvoi] = { 0.0 };
		eps_1[i] += D_EPS_CTAN_AVE;

		newton_raphson(false, eps_1, NULL, u_aux, NULL, NULL, NULL);

		calc_ave_stress(u_aux, NULL, sig_1);

		for (int v = 0; v < nvoi; ++v)
			ctan_lin[v * nvoi + i] = sig_1[v] / D_EPS_CTAN_AVE;
	}
	double max = ctan_lin[0];
	for (int i = 1; i < nvoi * nvoi; ++i)
		if (ctan_lin[i] > max)
			max = ctan_lin[i];

	for (int i = 0; i < nvoi * nvoi; ++i)
		ctan_lin[i] = (fabs(ctan_lin[i]) > max * 1.0e-3) ? ctan_lin[i] : 0.0;
}


template <int tdim>
material_t micropp<tdim>::get_material(const int e) const
{
	return material_list[elem_type[e]];
}


template <int tdim>
void micropp<tdim>::get_elem_rhs(double *u,
				 double *int_vars_old,
				 double be[npe * dim],
				 int ex, int ey, int ez) const
{
	constexpr int npedim = npe * dim;
	double bmat[nvoi][npedim], stress_gp[nvoi], strain_gp[nvoi];

	memset(be, 0, npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		calc_bmat(gp, bmat);

		get_strain(u, gp, strain_gp, ex, ey, ez);
		get_stress(gp, strain_gp, int_vars_old, stress_gp, ex, ey, ez);

		for (int i = 0; i < npedim; ++i)
			for (int j = 0; j < nvoi; ++j)
				be[i] += bmat[j][i] * stress_gp[j] * wg;
	}
}


template <int tdim>
void micropp<tdim>::get_elem_nodes(int n[npe], int ex, int ey, int ez) const
{
	const int nxny = ny * nx;
	const int n0 = ez * nxny + ey * nx + ex;
	n[0] = n0;
	n[1] = n0 + 1;
	n[2] = n0 + nx + 1;
	n[3] = n0 + nx;

	if (dim == 3) {
		n[4] = n[0] + nxny;
		n[5] = n[1] + nxny;
		n[6] = n[2] + nxny;
		n[7] = n[3] + nxny;
	}
}


template<int tdim>
int micropp<tdim>::get_elem_type(int ex, int ey, int ez) const
{
	assert(0 <= micro_type  && micro_type <=4);

	const double coor[3] = {
		ex * dx + dx / 2.,
		ey * dy + dy / 2.,
		ez * dz + dz / 2. }; // 2D -> dz = 0

	if (micro_type == 0) { // sphere in the center

		const double rad = special_param;
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		double tmp = 0.;
		for (int i = 0; i < dim; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == 1) { // 2 flat layers in y dir

		const double width = special_param;
		return (coor[1] < width);

	} else if (micro_type == 2) { // a cilindrical fiber in z dir

		const double rad = special_param;
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		double tmp = 0.;
		for (int i = 0; i < 2; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == 3) { // 2 cilindrical fibers one in x and z dirs

		const double rad = special_param;
		const double cen_1[3] = { lx / 2., ly * .75, lz / 2. };
		double tmp_1 = 0.;
		for (int i = 0; i < 2; ++i)
			tmp_1 += (cen_1[i] - coor[i]) * (cen_1[i] - coor[i]);

		const double cen_2[3] = { lx / 2., ly * .25, lz / 2. };
		double tmp_2 = 0.;
		for (int i = 1; i < 3; ++i)
			tmp_2 += (cen_2[i] - coor[i]) * (cen_2[i] - coor[i]);

		return ((tmp_1 < rad * rad) || (tmp_2 < rad * rad));

	} else if (micro_type == 4) { // 3 quad fibers in x, y and z dirs

		const double width = special_param;
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		return 0;
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}


template <int tdim>
void micropp<tdim>::get_elem_displ(const double *u,
				   double elem_disp[npe * dim],
				   int ex, int ey, int ez) const
{
	int n[npe] ;
	get_elem_nodes(n, ex, ey, ez);

	for (int i = 0 ; i < npe; ++i)
		for (int d = 0; d < dim; ++d)
			elem_disp[i * dim + d] = u[n[i] * dim + d];
}


template <int tdim>
void micropp<tdim>::get_strain(const double *u, int gp, double *strain_gp,
			       int ex, int ey, int ez) const
{
	double elem_disp[npe * dim];
	get_elem_displ(u, elem_disp, ex, ey, ez);

	double bmat[nvoi][npe * dim];
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
	printf("l = [%lf, %lf, %lf]; param = %lf\n", lx, ly, lz, special_param);
	for (int i = 0; i < numMaterials; ++i) {
		printf("Type = %d, E = %e, nu = %e, Sy = %e, Ka = %e, \
		       plast = %d\n",
		       material_list[i].type,
		       material_list[i].E, material_list[i].nu,
		       material_list[i].Sy, material_list[i].Ka,
		       material_list[i].plasticity);
	}
}


template <int tdim>
void micropp<tdim>::get_stress(int gp, const double eps[nvoi],
			       double *int_vars_old,
			       double stress_gp[nvoi],
			       int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t material = get_material(e);
	const double mu = material.mu;
	double zero_nvoi[nvoi] = { 0.0 };

	if (material.plasticity == true) {

		double *eps_p_old;
		double alpha_old ;

		if (int_vars_old != NULL) {
			eps_p_old = &int_vars_old[intvar_ix(e, gp, 0)];
			alpha_old = int_vars_old[intvar_ix(e, gp, 6)];
		} else {
			eps_p_old = zero_nvoi;
			alpha_old = 0.0;
		}

		plastic_get_stress(&material, eps, eps_p_old, alpha_old, stress_gp);

	} else {

		isolin_get_stress(&material, eps, stress_gp);
	}

}


template <int tdim>
void micropp<tdim>::calc_ave_stress(double *u, double *int_vars_old,
				    double stress_ave[nvoi]) const
{
	memset(stress_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double stress_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, int_vars_old,
						   stress_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						stress_aux[v] += stress_gp[v] * wg;

				}
				for (int v = 0; v < nvoi; ++v)
					stress_ave[v] += stress_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; ++v)
		stress_ave[v] /= vol_tot;
}


template <int tdim>
void micropp<tdim>::calc_ave_strain(const double *u,
				    double strain_ave[nvoi]) const
{
	memset(strain_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double strain_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {
					double strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						strain_aux[v] += strain_gp[v]\
								 * wg;
				}

				for (int v = 0; v < nvoi; v++)
					strain_ave[v] += strain_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] /= vol_tot;
}


template<int tdim>
void micropp<tdim>::calc_fields(double *u, double *int_vars_old)
{
	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double eps_a[nvoi] = { 0.0 };
				double sig_a[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, int_vars_old,
						   stress_gp, ex, ey, ez);

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


template<int tdim>
bool micropp<tdim>::calc_vars_new(double *u, double *int_vars_old,
				  double *int_vars_new)
{
	INST_START;

	bool nl_flag = false;
	double zero_nvoi[nvoi] = { 0.0 };

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex){

				const int e = glo_elem(ex, ey, ez);
				const material_t material = get_material(e);

				for (int gp = 0; gp < npe; ++gp) {

					if (material.plasticity == true) {

						double *eps_p_old;
						double alpha_old;
						if (int_vars_old != NULL) {
							eps_p_old = &int_vars_old[intvar_ix(e, gp, 0)];
							alpha_old = int_vars_old[intvar_ix(e, gp, 6)];
						} else {
							eps_p_old = zero_nvoi;
							alpha_old = 0.0;
						}
						double *eps_p_new = &int_vars_new[intvar_ix(e, gp, 0)];
						double *alpha_new = &int_vars_new[intvar_ix(e, gp, 6)];
						double eps[nvoi];
						get_strain(u, gp, eps, ex, ey, ez);

						nl_flag |= plastic_evolute(&material,
									   eps, eps_p_old,
									   alpha_old,
									   eps_p_new,
									   alpha_new);
					}
				}
			}
		}
	}

	return nl_flag;
}


template class micropp<2>;
template class micropp<3>;

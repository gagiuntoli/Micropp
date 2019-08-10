/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
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
#include "material.hpp"


template<int tdim>
micropp<tdim>::micropp(const micropp_params_t &params):

	ngp(params.ngp),
	nx(params.size[0]), ny(params.size[1]),
	nz((tdim == 3) ? params.size[2] : 1),

	nn(nx * ny * nz),
	nndim(nn * dim),

	nex(nx - 1), ney(ny - 1),
	nez((tdim == 3) ? (nz - 1) : 1),

	nelem(nex * ney * nez),
	lx(1.0), ly(1.0), lz((tdim == 3) ? 1.0 : 0.0),
	dx(lx / nex), dy(ly / ney), dz((tdim == 3) ? lz / nez : 0.0),

	subiterations(params.subiterations),
	nsubiterations(params.nsubiterations),
	mpi_rank(params.mpi_rank),

	wg(((tdim == 3) ? dx * dy * dz : dx * dy) / npe),
	vol_tot((tdim == 3) ? lx * ly * lz : lx * ly),
	ivol(1.0 / (wg * npe)),
	evol((tdim == 3) ? dx * dy * dz : dx * dy),
	micro_type(params.type), nvars(nelem * npe * NUM_VAR_GP),

	nr_max_its(params.nr_max_its),
	nr_max_tol(params.nr_max_tol),
	nr_rel_tol(params.nr_rel_tol),
	calc_ctan_lin_flag(params.calc_ctan_lin),

	use_A0(params.use_A0),
	its_with_A0(params.its_with_A0),
	lin_stress(params.lin_stress),
	write_log_flag(params.write_log)
{
	INST_CONSTRUCT; // Initialize the Intrumentation

	/* GPU device selection if they are accessible */

#ifdef _OPENACC
	int acc_num_gpus = acc_get_num_devices(acc_device_nvidia);
	gpu_id = mpi_rank % acc_num_gpus;
	acc_set_device_num(gpu_id, acc_device_nvidia);
#endif

	for (int gp = 0; gp < npe; gp++) {
		calc_bmat(gp, calc_bmat_cache[gp]);
	}

	gp_list = new gp_t<tdim>[ngp]();
	for (int gp = 0; gp < ngp; ++gp) {
		if(params.coupling != nullptr) {
			gp_list[gp].coupling = params.coupling[gp];
			switch (params.coupling[gp]) {
				case FE_LINEAR:
					num_fe_linear ++;
					break;
				case FE_ONE_WAY:
					num_fe_one_way ++;
					break;
				case FE_FULL:
					num_fe_full ++;
					break;
			}
		} else {
			gp_list[gp].coupling = FE_ONE_WAY;
			num_fe_one_way ++;
		}

		num_fe_points = num_fe_linear + num_fe_one_way + num_fe_full;

		gp_list[gp].nndim = nndim;
		gp_list[gp].nvars = nvars;

		if (params.coupling == nullptr ||
		    (params.coupling[gp] == FE_ONE_WAY || params.coupling[gp] == FE_FULL)) {
			gp_list[gp].allocate_u();
		}
	}

	elem_type = (int *) calloc(nelem, sizeof(int));
	elem_stress = (double *) calloc(nelem * nvoi, sizeof(double));
	elem_strain = (double *) calloc(nelem * nvoi, sizeof(double));

	for (int i = 0; i < num_geo_params; i++) {
		geo_params[i] = params.geo_params[i];
	}

	for (int i = 0; i < MAX_MATERIALS; ++i) {
		material_list[i] = material_t::make_material(params.materials[i]);
	}

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {
				const int e_i = glo_elem(ex, ey, ez);
				elem_type[e_i] = get_elem_type(ex, ey, ez);
			}
		}
	}

	calc_volume_fractions();

	if (params.use_A0) {
#ifdef _OPENMP
		int num_of_A0s = omp_get_max_threads();
#else
		int num_of_A0s = 1;
#endif
		A0 = (ell_matrix *) malloc(num_of_A0s * sizeof(ell_matrix));

#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < num_of_A0s; ++i) {
			ell_init(&A0[i], dim, dim, params.size, CG_ABS_TOL, CG_REL_TOL, CG_MAX_ITS);
			double *u = (double *) calloc(nndim, sizeof(double));
			assembly_mat(&A0[i], u, nullptr);
			free(u);
		}
	}

	/* Average tangent constitutive tensor initialization */

	memset(ctan_lin_fe, 0.0, nvoi * nvoi * sizeof(double));

	if (calc_ctan_lin_flag) {
		if (num_fe_points > 0) {
			calc_ctan_lin_fe_models();
		}
	}

	for (int gp = 0; gp < ngp; ++gp) {

		if (gp_list[gp].coupling == FE_LINEAR || gp_list[gp].coupling == FE_ONE_WAY ||
		    gp_list[gp].coupling == FE_FULL) {

			memcpy(gp_list[gp].ctan, ctan_lin_fe, nvoi * nvoi * sizeof(double));

		} else if (gp_list[gp].coupling == MIXTURE_RULE_CHAMIS) {

			double ctan[nvoi * nvoi];
			calc_ctan_lin_mixture_rule_Chamis(ctan);
			memcpy(gp_list[gp].ctan, ctan, nvoi * nvoi * sizeof(double));

		}
	}

	/* Open the log file */

	if (write_log_flag) {
		char filename[128];
		std::stringstream filename_stream;
		filename_stream << "micropp-profiling-" << mpi_rank << ".log";
		std::string file_name_string = filename_stream.str();
		strcpy(filename, file_name_string.c_str());

		ofstream_log.open(filename, ios::out);
		ofstream_log << "#<gp_id>  <non-linear>  <cost>  <converged>" << endl;
	}

}


template <int tdim>
micropp<tdim>::~micropp()
{
	INST_DESTRUCT;

	free(elem_stress);
	free(elem_strain);
	free(elem_type);

	if (use_A0) {
#ifdef _OPENMP
		int num_of_A0s = omp_get_max_threads();
#else
		int num_of_A0s = 1;
#endif

#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < num_of_A0s; ++i) {
			ell_free(&A0[i]);
		}
		free(A0);
	}

	for (int i = 0; i < MAX_MATERIALS; ++i) {
		delete material_list[i];
	}

	delete [] gp_list;
}


template <int tdim>
int micropp<tdim>::is_non_linear(const int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return (int) gp_list[gp_id].allocated;
}


template <int tdim>
int micropp<tdim>::get_cost(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].cost;
}


template <int tdim>
bool micropp<tdim>::has_converged(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].converged;
}


template <int tdim>
bool micropp<tdim>::has_subiterated(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].subiterated;
}


template <int tdim>
int micropp<tdim>::get_non_linear_gps(void) const
{
	int count = 0;
	for (int gp = 0; gp < ngp; ++gp) {
		if (gp_list[gp].allocated) {
			count ++;
		}
	}
	return count;
}


template <int tdim>
void micropp<tdim>::calc_ctan_lin_fe_models()
{

#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < nvoi; ++i) {

		const int ns[3] = { nx, ny, nz };

		ell_matrix A;  // Jacobian
		ell_init(&A, dim, dim, ns, CG_ABS_TOL, CG_REL_TOL, CG_MAX_ITS);
		double *b = (double *) calloc(nndim, sizeof(double));
		double *du = (double *) calloc(nndim, sizeof(double));
		double *u = (double *) calloc(nndim, sizeof(double));

		double sig[6];
		double eps[nvoi] = { 0.0 };
		eps[i] += D_EPS_CTAN_AVE;

		newton_raphson(&A, b, u, du, eps);

		calc_ave_stress(u, sig);

		for (int v = 0; v < nvoi; ++v) {
			ctan_lin_fe[v * nvoi + i] = sig[v] / D_EPS_CTAN_AVE;
		}

		ell_free(&A);
		free(b);
		free(u);
		free(du);
	}
}


template <int tdim>
void micropp<tdim>::calc_ctan_lin_mixture_rule_Chamis(double ctan[nvoi * nvoi])
{

 const double Em = material_list[0]->E;
 const double nu_m = material_list[0]->nu;
 const double Gm = Em / (2 * (1 + nu_m));

 const double Ef = material_list[1]->E;
 const double nu_f = material_list[1]->nu;
 const double Gf = Ef / (2 * (1 + nu_f));

 const double E11 = Vf * Ef + Vm * Em;
 const double E22 = Em / (1 - sqrt(Vf) * (1 - Em / Ef));
 const double nu_12 = Vf * nu_f + Vm * nu_m;
 const double G12 = Gm / (1 - sqrt(Vf) * (1 - Gm / Gf));

 const double S[3][3] = {
	 {       1 / E11, - nu_12 / E11, - nu_12 / E11 },
	 { - nu_12 / E11,       1 / E22, - nu_12 / E22 },
	 { - nu_12 / E11, - nu_12 / E22,       1 / E22 },
 };

 double det =
	 S[0][0] * (S[1][1] * S[2][2] - S[2][1] * S[1][2]) -
	 S[0][1] * (S[1][0] * S[2][2] - S[2][0] * S[1][2]) +
	 S[0][2] * (S[1][0] * S[2][1] - S[2][0] * S[1][1]);

 double c00 = +(S[1][1] * S[2][2] - S[2][1] * S[1][2]);
 double c01 = -(S[1][0] * S[2][2] - S[2][0] * S[1][2]);
 double c02 = +(S[1][0] * S[2][1] - S[2][0] * S[1][1]);

 double c10 = -(S[0][1] * S[2][2] - S[2][1] * S[0][2]);
 double c11 = +(S[0][0] * S[2][2] - S[2][0] * S[0][2]);
 double c12 = -(S[0][0] * S[2][1] - S[2][0] * S[0][1]);

 double c20 = +(S[0][1] * S[1][2] - S[1][1] * S[0][2]);
 double c21 = -(S[0][0] * S[1][2] - S[1][0] * S[0][2]);
 double c22 = +(S[0][0] * S[1][1] - S[1][0] * S[0][1]);

 memset (ctan, 0, nvoi * nvoi * sizeof(double));

 ctan[0 * nvoi + 0] = c00 / det;
 ctan[0 * nvoi + 1] = c01 / det;
 ctan[0 * nvoi + 2] = c02 / det;

 ctan[1 * nvoi + 0] = c10 / det;
 ctan[1 * nvoi + 1] = c11 / det;
 ctan[1 * nvoi + 2] = c12 / det;

 ctan[2 * nvoi + 0] = c20 / det;
 ctan[2 * nvoi + 1] = c21 / det;
 ctan[2 * nvoi + 2] = c22 / det;

 ctan[3 * nvoi + 3] = G12;
 ctan[4 * nvoi + 4] = G12;
 ctan[5 * nvoi + 5] = G12;

}


template <int tdim>
material_t *micropp<tdim>::get_material(const int e) const
{
	return material_list[elem_type[e]];
}


template <int tdim>
void micropp<tdim>::get_elem_rhs(const double *u,
				 const double *vars_old,
				 double be[npe * dim],
				 int ex, int ey, int ez) const
{
	constexpr int npedim = npe * dim;
	double stress_gp[nvoi], strain_gp[nvoi];

	memset(be, 0, npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		get_strain(u, gp, strain_gp, ex, ey, ez);
		get_stress(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);

		for (int i = 0; i < npedim; ++i)
			for (int j = 0; j < nvoi; ++j)
				be[i] += calc_bmat_cache[gp][j][i] * stress_gp[j] * wg;
	}
}


template <int tdim>
void micropp<tdim>::get_elem_mat(const double *u,
				 const double *vars_old,
				 double Ae[npe * dim * npe * dim],
				 int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t *material = get_material(e);

	double ctan[nvoi][nvoi];
	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;

	double TAe[npedim2] = { 0.0 };

	for (int gp = 0; gp < npe; ++gp) {

		double eps[6];
		get_strain(u, gp, eps, ex, ey, ez);

		const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;
		material->get_ctan(eps, (double *)ctan, vars);

		double cxb[nvoi][npedim];

		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < npedim; ++j) {
				double tmp = 0.0;
				for (int k = 0; k < nvoi; ++k)
					tmp += ctan[i][k] * calc_bmat_cache[gp][k][j];
				cxb[i][j] = tmp * wg;
			}
		}

		for (int m = 0; m < nvoi; ++m) {
			for (int i = 0; i < npedim; ++i) {
				const int inpedim = i * npedim;
				const double bmatmi = calc_bmat_cache[gp][m][i];
				for (int j = 0; j < npedim; ++j)
					TAe[inpedim + j] += bmatmi * cxb[m][j];
			}
		}
	}
	memcpy(Ae, TAe, npedim2 * sizeof(double));
}


#pragma acc routine seq
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
	const double coor[3] = {
		ex * dx + dx / 2.,
		ey * dy + dy / 2.,
		ez * dz + dz / 2. }; // 2D -> dz = 0

	if (micro_type == MIC_HOMOGENEOUS) { // Only one material (mat[0])

		return 0;

	} else if (micro_type == MIC_SPHERE) { // sphere in the center

		const double rad = geo_params[0];
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0

		return point_inside_sphere(center, rad, coor);

	} else if (micro_type == MIC_LAYER_Y) { // 2 flat layers in y dir

		const double width = geo_params[0];
		return (coor[1] < width);

	} else if (micro_type == MIC_CILI_FIB_X) { // a cilindrical fiber in x dir

		const double rad = geo_params[0];
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		const double dir[3] = { 1, 0, 0 };

		return point_inside_cilinder_inf(dir, center, rad, coor);

	} else if (micro_type == MIC_CILI_FIB_Z) { // a cilindrical fiber in z dir

		const double rad = geo_params[0];
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		const double dir[3] = { 0, 0, 1 };

		return point_inside_cilinder_inf(dir, center, rad, coor);

	} else if (micro_type == MIC_CILI_FIB_XZ) { // 2 cilindrical fibers one in x and z dirs

		const double rad = geo_params[0];
		const double cen_1[3] = { lx / 2., ly * .75, lz / 2. };
		const double cen_2[3] = { lx / 2., ly * .25, lz / 2. };
		const double dir_x[3] = { 1, 0, 0 };
		const double dir_z[3] = { 0, 0, 1 };

		if (point_inside_cilinder_inf(dir_z, cen_1, rad, coor) ||
		    point_inside_cilinder_inf(dir_x, cen_2, rad, coor)) {
			return 1;
		}
		return 0;

	} else if (micro_type == MIC_QUAD_FIB_XYZ) {

	       	/* 3 quad fibers in x, y and z dirs */

		const double width = geo_params[0];
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

	} else if (micro_type == MIC_QUAD_FIB_XZ) {

	       	/* 2 quad fibers in x and z dirs */

		const double width = geo_params[0];
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		return 0;

	} else if (micro_type == MIC_QUAD_FIB_XZ_BROKEN_X) {

	       	/* 2 quad fibers in x and z dirs and the one in x is broken */

		const double width = geo_params[0];
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width &&
		    (coor[0] < lx * .8 || coor[0] > lx * .9))
			return 1;

		return 0;

	} else if (micro_type == MIC3D_SPHERES) {

	       	/* Distribution of Several Spheres of diferent sizes */

		const int num_spheres = 40;

		const double centers[num_spheres][3] = {
			{ .8663, .0689, .1568 },
			{ .2305, .4008, .2093 },
			{ .1987, .8423, .2126 },
			{ .6465, .7095, .4446 },
			{ .6151, .9673, .2257 },
			{ .1311, .9739, .4129 },
			{ .8433, .0738, .2233 },
			{ .5124, .2111, .0369 },
			{ .4094, .7030, .6241 },
			{ .1215, .8289, .3812 },
			{ .1125, .9266, .6872 },
			{ .7422, .4030, .6000 },
			{ .0791, .3819, .9297 },
			{ .1201, .7712, .0069 },
			{ .9339, .2177, .1976 },
			{ .3880, .1331, .9032 },
			{ .7319, .7146, .8150 },
			{ .8796, .2858, .6701 },
			{ .1427, .9309, .9830 },
			{ .1388, .3126, .8054 },
			{ .0729, .4754, .4950 },
			{ .5038, .7755, .7333 },
			{ .0242, .2185, .9884 },
			{ .5371, .5372, .8269 },
			{ .4564, .0039, .1715 },
			{ .7410, .0799, .8775 },
			{ .2627, .9047, .5681 },
			{ .7894, .1908, .2993 },
			{ .1373, .5370, .9916 },
			{ .0207, .8020, .9283 },
			{ .6540, .2359, .0286 },
			{ .7801, .3568, .5086 },
			{ .7322, .5706, .1114 },
			{ .5479, .1493, .1267 },
			{ .6722, .1530, .1003 },
			{ .7659, .3426, .9181 },
			{ .4582, .4636, .6310 },
			{ .1811, .9665, .4713 },
			{ .0834, .6066, .6936 },
			{ .4865, .7584, .3635 }
		};

		const double rads[num_spheres] = {
			0.1 * .5741,
			0.1 * .1735,
			0.1 * .5065,
			0.1 * .8565,
			0.1 * .9735,
			0.1 * .7087,
			0.1 * .9585,
			0.1 * .2843,
			0.1 * .4029,
			0.1 * .2574,
			0.1 * .1575,
			0.1 * .3316,
			0.1 * .1874,
			0.1 * .6300,
			0.1 * .7049,
			0.1 * .7258,
			0.1 * .8002,
			0.1 * .8176,
			0.1 * .9970,
			0.1 * .4736,
			0.1 * .6351,
			0.1 * .3140,
			0.1 * .3440,
			0.1 * .2444,
			0.1 * .3894,
			0.1 * .5234,
			0.1 * .1658,
			0.1 * .4020,
			0.1 * .6813,
			0.1 * .2060,
			0.1 * .1427,
			0.1 * .9332,
			0.1 * .5406,
			0.1 * .9765,
			0.1 * .0956,
			0.1 * .6432,
			0.1 * .9998,
			0.1 * .4166,
			0.1 * .6907,
			0.1 * .0404
		};


		for (int i = 0; i < num_spheres; ++i) {
			if (point_inside_sphere(centers[i], rads[i], coor)) {
				return 1;
			}
		}

		return 0;

	} else if (micro_type == MIC3D_8) {

		/*
		 * returns
		 * 0 : for the matrix
		 * 1 : for the cilinders
		 * 2 : for the layer around the cilinders
		 * 2 : for the flat layer
		 *
		 */

		const double rad_cilinder = geo_params[0];
		const double width_flat_layer = geo_params[1];
		const double width_cili_layer = geo_params[2];

		const double cen_1[3] = { lx * .25, ly * .75, 0.0 };
		const double cen_2[3] = { lx * .75, ly * .75, 0.0 };

		const double cen_3[3] = { 0.0, ly * .25, lz * .25 };
		const double cen_4[3] = { 0.0, ly * .25, lz * .75 };

		const double dir_x[3] = { 1, 0, 0 };
		const double dir_z[3] = { 0, 0, 1 };

		if(point_inside_cilinder_inf(dir_z, cen_1, rad_cilinder, coor) ||
		   point_inside_cilinder_inf(dir_z, cen_2, rad_cilinder, coor) ||
		   point_inside_cilinder_inf(dir_x, cen_3, rad_cilinder, coor) ||
		   point_inside_cilinder_inf(dir_x, cen_4, rad_cilinder, coor)) {
			return 1;
		}

		if(point_inside_cilinder_inf(dir_z, cen_1, rad_cilinder + width_cili_layer, coor) ||
		   point_inside_cilinder_inf(dir_z, cen_2, rad_cilinder + width_cili_layer, coor) ||
		   point_inside_cilinder_inf(dir_x, cen_3, rad_cilinder + width_cili_layer, coor) ||
		   point_inside_cilinder_inf(dir_x, cen_4, rad_cilinder + width_cili_layer, coor) ||
		   fabs(coor[1] - ly / 2) < width_flat_layer) {
			return 2;
		}

		return 0;


	} else if (micro_type == MIC3D_FIBS_20_DISORDER) {

		const int num_fibs = 20;
		const double radius = 0.05;

		const double dirs[num_fibs][3] = {
			{ 1, .5741, .8515 },
			{ 1, .1735, .1103 },
			{ 1, .5065, .4600 },
			{ 1, .8565, .9045 },
			{ 1, .9735, .6313 },
			{ 1, .7087, .1547 },
			{ 1, .9585, .0220 },
			{ 1, .2843, .4062 },
			{ 1, .4029, .8095 },
			{ 1, .2574, .4742 },
			{ 1, .1575, .0768 },
			{ 1, .3316, .0320 },
			{ 1, .1874, .6364 },
			{ 1, .6300, .2688 },
			{ 1, .7049, .5137 },
			{ 1, .7258, .2799 },
			{ 1, .8002, .4794 },
			{ 1, .8176, .4142 },
			{ 1, .9970, .2189 },
			{ 1, .4736, .6202 }
		};


		const double centers[num_fibs][3] = {
			{ .8663, .0689, .1568 },
			{ .2305, .4008, .2093 },
			{ .1987, .8423, .2126 },
			{ .6465, .7095, .4446 },
			{ .6151, .9673, .2257 },
			{ .1311, .9739, .4129 },
			{ .8433, .0738, .2233 },
			{ .5124, .2111, .0369 },
			{ .4094, .7030, .6241 },
			{ .1215, .8289, .3812 },
			{ .1125, .9266, .6872 },
			{ .7422, .4030, .6000 },
			{ .0791, .3819, .9297 },
			{ .1201, .7712, .0069 },
			{ .9339, .2177, .1976 },
			{ .3880, .1331, .9032 },
			{ .7319, .7146, .8150 },
			{ .8796, .2858, .6701 },
			{ .1427, .9309, .9830 },
			{ .1388, .3126, .8054 }
		};

		for (int i = 0; i < num_fibs; ++i) {
			if(point_inside_cilinder_inf(dirs[i], centers[i], radius, coor)) {
				return 1;
			}
		}
		return 0;

	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}


#pragma acc routine seq
template <int tdim>
void micropp<tdim>::get_elem_displ(const double *u,
				   double elem_disp[npe * dim],
				   int ex, int ey, int ez) const
{
	int n[npe] ;
	get_elem_nodes(n, ex, ey, ez);

	for (int i = 0 ; i < npe; ++i) {
		for (int d = 0; d < dim; ++d) {
			elem_disp[i * dim + d] = u[n[i] * dim + d];
		}
	}
}


#pragma acc routine seq
template <int tdim>
void micropp<tdim>::get_strain(const double *u, int gp, double *strain_gp,
			       int ex, int ey, int ez) const
{
	double elem_disp[npe * dim];
	get_elem_displ(u, elem_disp, ex, ey, ez);

	for (int i = 0; i < nvoi; ++i) {
		strain_gp[i] = 0;
	}

	for (int v = 0; v < nvoi; ++v) {
		for (int i = 0; i < npe * dim; i++){
			strain_gp[v] += calc_bmat_cache[gp][v][i] * elem_disp[i];
		}
	}
}


template <int tdim>
void micropp<tdim>::print_info() const
{
	cout << "micropp" << dim << endl;

	cout << "Micro-structure   : " << micro_names[micro_type] << endl;

	cout << "MATRIX [%]        : " << Vm << endl;
	cout << "FIBER  [%]        : " << Vf << endl;

	cout << "FE_LINEAR         : " << num_fe_linear    << " GPs" << endl;
	cout << "FE_ONE_WAY        : " << num_fe_one_way   << " GPs" << endl;
	cout << "FE_FULL           : " << num_fe_full      << " GPs" << endl;
	cout << "USE A0            : " << use_A0 << endl;
	cout << "NUM SUBITS        : " << nsubiterations << endl;
	cout << "MPI RANK          : " << mpi_rank << endl;

#ifdef _OPENACC
	int acc_num_gpus = acc_get_num_devices(acc_device_nvidia);
	cout << "ACC NUM GPUS      : " << acc_num_gpus << endl;
	cout << "GPU ID            : " << gpu_id << endl;
#endif
#ifdef _OPENMP
	int omp_max_threads = omp_get_max_threads();
	cout << "OMP NUM THREADS   : " << omp_max_threads << endl;
#endif

	cout    << "ngp :" << ngp
		<< " nx :" << nx << " ny :" << ny << " nz :" << nz
		<< " nn :" << nn << endl
		<< "lx : " << lx << " ly : " << ly << " lz : " << lz;
	cout << endl;

	cout << "geo_params:";
	for (int i = 0; i < num_geo_params; ++i) {
		cout << " " << geo_params[i];
	}
	cout << endl;

	for (int i = 0; i < MAX_MATERIALS; ++i) {
		material_list[i]->print();
		cout << endl;
	}

	cout << endl;
	cout << "ctan_lin_fe = " << endl;
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			cout << ctan_lin_fe[i * 6 + j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}


template <int tdim>
void micropp<tdim>::get_stress(int gp, const double eps[nvoi],
			       const double *vars_old,
			       double stress_gp[nvoi],
			       int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t *material = get_material(e);
	const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;

	material->get_stress(eps, stress_gp, vars);
}


template <int tdim>
void micropp<tdim>::calc_ave_stress(const double *u, double stress_ave[nvoi],
				    const double *vars_old) const
{
	memset(stress_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double stress_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];
					get_strain(u, gp, strain_gp, ex, ey, ez);
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

	for (int v = 0; v < nvoi; ++v)
		stress_ave[v] /= vol_tot;
}


template <int tdim>
void micropp<tdim>::calc_volume_fractions()
{
	Vm = 0.0;
	Vf = 0.0;
	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {
				const int e_i = glo_elem(ex, ey, ez);
				if (elem_type[e_i] == 0) {
					Vm += evol;
				} else if (elem_type[e_i] >= 1) {
					Vf += evol;
				}
			}
		}
	}
	Vm /= vol_tot;
	Vf /= vol_tot;
}


template <int tdim>
void micropp<tdim>::calc_ave_strain(const double *u, double strain_ave[nvoi]) const
{
	memset(strain_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double strain_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {
					double strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
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


template<int tdim>
void micropp<tdim>::calc_fields(double *u, double *vars_old)
{
	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double eps_a[nvoi] = { 0.0 };
				double sig_a[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
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
					get_strain(u, gp, eps, ex, ey, ez);

					non_linear |= material->evolute(eps, vars_old, vars_new);
				}
			}
		}
	}

	return non_linear;
}


template class micropp<3>;

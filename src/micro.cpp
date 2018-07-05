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

#include <vector>
#include <iostream>

#include <cassert>

#include "micro.hpp"

micropp_t::micropp_t(const int _dim, const int size[3], const int _micro_type,
                     const double *_micro_params, const int *mat_types, const double *params):
	dim(_dim),

	lx(_micro_params[0]),
	ly(_micro_params[1]),
	lz(dim == 2 ? 0.0 : _micro_params[2]),
	width(_micro_params[3]),

	inv_tol(micro_params[4]),
	nx(size[0]),
	ny(size[1]),
	nz(_dim == 2 ? 1 : size[2]),
	nn(nx * ny * nz),

	dx(lx / (nx - 1)),
	dy(ly / (ny - 1)),
	dz(lz / (nz - 1)),

	npe(dim == 2 ? 4 : 8),
	nvoi(dim == 2 ? 3 : 6),
	nelem(dim == 2 ? (nx - 1) * (ny - 1) : (nx - 1) * (ny - 1) * (nz - 1)),

	micro_type(_micro_type),
	num_int_vars(nelem * 8 * NUM_VAR_GP)
{
	assert(dim == 2 || dim == 3);


	b = (double *) malloc(nn * dim * sizeof(double));
	du = (double *) malloc(nn * dim * sizeof(double));
	u = (double *) malloc(nn * dim * sizeof(double));
	elem_stress = (double *) malloc(nelem * nvoi * sizeof(double));
	elem_strain = (double *) malloc(nelem * nvoi * sizeof(double));
	elem_type = (int *) malloc(nelem * sizeof(int));
	vars_old = (double *) malloc(num_int_vars * sizeof(double));
	vars_new = (double *) malloc(num_int_vars * sizeof(double));

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
		material_list[i].E  = params[i * MAX_MAT_PARAM + 0];
		material_list[i].nu = params[i * MAX_MAT_PARAM + 1];
		material_list[i].Sy = params[i * MAX_MAT_PARAM + 2];
		material_list[i].Ka = params[i * MAX_MAT_PARAM + 3];

		material_list[i].k = material_list[i].E / (3 * (1 - 2 * material_list[i].nu));
		material_list[i].mu = material_list[i].E / (2 * (1 + material_list[i].nu));
		material_list[i].lambda = (material_list[i].nu * material_list[i].E) /
			((1 + material_list[i].nu) * (1 - 2 * material_list[i].nu));	// lambda

		if (mat_types[i] == 0) {
			// lineal
			material_list[i].plasticity = false;
			material_list[i].damage = false;
		} else if (mat_types[i] == 1) {
			// con plasticidad
			material_list[i].plasticity = true;
			material_list[i].damage = false;
		} else if (mat_types[i] == 2) {
			// con daño
			material_list[i].plasticity = false;
			material_list[i].damage = true;
		}
		file << setw(14) << material_list[i].E << " " << material_list[i].nu
		     << " " << material_list[i].Sy << " " << material_list[i].Ka << endl;
	}
	file.close();

	if (dim == 2) {
		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				int e = glo_elem3D(ex, ey, 0);
				elem_type[e] = get_elem_type2D(ex, ey);
			}
		}
		ell_init_2D(&A, dim, nx, ny);
	} else if (dim == 3) {
		for (int ex = 0; ex < nx - 1; ex++) {
			for (int ey = 0; ey < ny - 1; ey++) {
				for (int ez = 0; ez < nz - 1; ez++) {
					int e = glo_elem3D(ex, ey, ez);
					elem_type[e] = get_elem_type3D(ex, ey, ez);
				}
			}
		}
		ell_init_3D(&A, dim, nx, ny, nz);
	}

	calc_ctan_lin();

	output_files_header = false;

	file.open("micropp_convergence.dat");
	file.close();
	file.open("micropp_eps_sig_ctan.dat");
	file.close();
	file.open("micropp_int_vars_n.dat");
	file.close();
}

micropp_t::~micropp_t()
{
	ell_free(&A);

	free(b);
	free(du);
	free(u);
	free(elem_stress);
	free(elem_strain);
	free(elem_type);
	free(vars_old);
	free(vars_new);

	for (auto const &gp:gauss_list) {
		free(gp.int_vars_n);
		free(gp.int_vars_k);
	}
}

void micropp_t::get_nl_flag(int gp_id, int *non_linear)
{
	*non_linear = 0;
	for (auto const &gp : gauss_list)
		if (gp.id == gp_id) {
			*non_linear = (gp.int_vars_n != NULL);
			break;
		}
}

int micropp_t::get_elem_type2D(int ex, int ey)
{
	assert(micro_type == 0 || micro_type == 1);

	if (micro_type == 0) {
		// esfera en matriz
		double x1 = ex * dx + dx / 2;
		double y1 = ey * dy + dy / 2;
		double x2 = lx / 2;
		double y2 = ly / 2;
		double rad = micro_params[3];
		return ((x2 - x1) * (x2 - x1) +
		        (y2 - y1) * (y2 - y1) < width * width);

	} else if (micro_type == 1) {
		const double y = ey * dy + dy / 2;
		return (y < width);
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}

int micropp_t::get_elem_type3D(int ex, int ey, int ez)
{
	assert(micro_type == 0 || micro_type == 1);

	if (micro_type == 0) {
		// esfera en matriz
		double x1 = ex * dx + dx / 2;
		double y1 = ey * dy + dy / 2;
		double z1 = ez * dz + dz / 2;
		double x2 = lx / 2;
		double y2 = ly / 2;
		double z2 = lz / 2;
		double rad = micro_params[3];
		return ((x2 - x1) * (x2 - x1) +
		        (y2 - y1) * (y2 - y1) +
		        (z2 - z1) * (z2 - z1) < width * width);

	} else if (micro_type == 1) {
		const double y = ey * dy + dy / 2;
		return (y < width);
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}

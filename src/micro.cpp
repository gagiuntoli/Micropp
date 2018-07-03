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
#include "micro.h"

Problem::Problem (int dim, int size[3], int micro_type, double *micro_params, int *mat_types, double *params)
{
  	this->micro_type = micro_type;

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

  	for (int i=0; i<nParams; i++) {
    	this->micro_params[i] = micro_params[i];
  	}
  	lx = this->micro_params[0];
  	ly = this->micro_params[1];
  	lz = this->micro_params[2];

  	inv_tol = this->micro_params[4];  
  	
  	ofstream file;
  	file.open ("micropp_materials.dat");
    file << scientific;
    	
  	for (int i=0; i<numMaterials; i++)
  	{
    	material_list[i].E  = params[i*MAX_MAT_PARAM + 0];
    	material_list[i].nu = params[i*MAX_MAT_PARAM + 1];
    	material_list[i].Sy = params[i*MAX_MAT_PARAM + 2];
    	material_list[i].Ka = params[i*MAX_MAT_PARAM + 3];

    	material_list[i].k  = material_list[i].E /(3*(1-2*material_list[i].nu));
    	material_list[i].mu = material_list[i].E /(2*(1+material_list[i].nu));
    	material_list[i].lambda = (material_list[i].nu*material_list[i].E) /((1+material_list[i].nu)*(1-2*material_list[i].nu)); // lambda

    	if (mat_types[i] == 0) {
      		// lineal
      		material_list[i].plasticity = false;
      		material_list[i].damage     = false;
    	} else if (mat_types[i] == 1) {
      		// con plasticidad
      		material_list[i].plasticity = true;
      		material_list[i].damage     = false;
    	} else if (mat_types[i] == 2) {
      		// con daño
      		material_list[i].plasticity = false;
      		material_list[i].damage     = true;
    	}
    	file << setw(14) << material_list[i].E << " " 
			 << material_list[i].nu << " " 
			 << material_list[i].Sy << " " 
			 << material_list[i].Ka << endl;
  	}
  	file.close ();

  	this->dim = dim;
  	if (dim == 2) {
    	nx = size[0];
    	ny = size[1];
    	nvoi = 3;
    	nn = nx * ny;
    	npe = 4;
  	} else if (dim == 3) {
    	nx = size[0];
    	ny = size[1];
    	nz = size[2];
    	nvoi = 6;
    	nn = nx * ny * nz;
    	npe = 8;
  	}

  	dx = lx/(nx-1);
  	dy = ly/(ny-1);
  	if (dim == 2) {
    	nelem = (nx-1) * (ny-1);
  	} else if (dim == 3) {
    	dz = lz/(nz-1);
    	nelem = (nx-1) * (ny-1) * (nz-1);
  	}

  	num_int_vars = nelem * 8 * NUM_VAR_GP;

  	b  = (double*)malloc(nn*dim*sizeof(double));
  	du = (double*)malloc(nn*dim*sizeof(double));
  	u  = (double*)malloc(nn*dim*sizeof(double));
  	elem_stress = (double*)malloc(nelem*nvoi*sizeof(double));
  	elem_strain = (double*)malloc(nelem*nvoi*sizeof(double));
  	elem_type = (int*)malloc(nelem*sizeof(int));
  	vars_old = (double*)malloc(num_int_vars*sizeof(double));
  	vars_new = (double*)malloc(num_int_vars*sizeof(double));

  	for (int i=0; i<nn*dim; i++)
    	u[i] = 0.0;

  	if (dim == 2) {
    	for (int ex=0; ex<nx-1; ex++) {
      		for (int ey=0; ey<ny-1; ey++) {
				int e = glo_elem3D(ex,ey,0);
				elem_type[e] = getElemType(ex,ey);
      		}
    	}
  	} else {
    	for (int ex=0; ex<nx-1; ex++) {
      		for (int ey=0; ey<ny-1; ey++) {
				for (int ez=0; ez<nz-1; ez++) {
	  				int e = glo_elem3D(ex,ey,ez);
	  				elem_type[e] = getElemType(ex,ey,ez);
				}
      		}
    	}
  	}

  	if (dim == 2)
    	ell_init_2D (A, dim, nx, ny);
  	else if (dim == 3)
    	ell_init_3D (A, dim, nx, ny, nz);

  	calcCtanLinear ();

  	output_files_header = false;

  	file.open ("micropp_convergence.dat");
  	file.close ();
  	file.open ("micropp_eps_sig_ctan.dat");
  	file.close ();
}

Problem::~Problem()
{
  	ell_free (A);
  	free(b);
  	free(du);
  	free(u);
  	free(elem_stress);
  	free(elem_strain);
  	free(elem_type);
  	free(vars_old);
  	free(vars_new);

	for (auto const& gp : gauss_list) {
    	free(gp.int_vars_n);
    	free(gp.int_vars_k);
  	}
}

void Problem::getNonLinearFlag (int gp_id, int *non_linear)
{
  	*non_linear = 0;
	for (auto const& gp : gauss_list)
    	if (gp.id == gp_id) {
      		*non_linear = (gp.int_vars_n == NULL) ? 0:1;
      		break;
    	}
}

int Problem::getElemType (int ex, int ey)
{
  	if (micro_type == 0) {
    	// esfera en matriz
    	double x1 = ex*dx + dx/2;
    	double y1 = ey*dy + dy/2;
    	double x2 = lx/2;
    	double y2 = ly/2;
    	double rad = micro_params[3];
    	if ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) < rad*rad)
      		return 1;
    	else 
      		return 0;

  	} else if (micro_type == 1) {
    	double y = ey*dy + dy/2;
    	double espesor = micro_params[3];
    	if (y < espesor)
      		return 1;
    	else 
      		return 0;
  	}
}

int Problem::getElemType (int ex, int ey, int ez)
{
  	if (micro_type == 0) {
    	// esfera en matriz
    	double x1 = ex*dx + dx/2;
    	double y1 = ey*dy + dy/2;
    	double z1 = ez*dz + dz/2;
    	double x2 = lx/2;
    	double y2 = ly/2;
    	double z2 = lz/2;
    	double rad = micro_params[3];
    	if ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) < rad*rad)
      		return 1;
    	else 
      		return 0;

  	} else if (micro_type == 1) {
    	double y = ey*dy + dy/2;
    	double espesor = micro_params[3];
    	if (y < espesor)
      		return 1;
    	else 
      		return 0;
  	}
}

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


#pragma once


#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>


using namespace std;


#define CG_ABS_TOL      1.0e-50
#define CG_MAX_ITS      1000
#define CG_REL_TOL      1.0e-5

#define nod_index(i,j,k)   ((k)*nx*ny + (j)*nx + (i))
#define nod_index3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define nod_index2D(i,j)   ((j)*nx + (i))


typedef struct {

	int n[3];          // nx ny nz
	int nn;
	int dim;
	int nfield;
	int shift;
	int nrow;          // number of rows
	int ncol;          // number of columns
	int nnz;           // non zeros per row
	int *cols = NULL;
	double *vals = NULL;

	int max_its;       // maximun number of iterations
	double min_err;    // minimun error (absolute)
	double rel_err;    // relative error
	double *k, *r, *z, *p, *Ap;

} ell_matrix;

void ell_init(ell_matrix *m, const int nfield, const int dim, const int ns[3],
	      const double min_err = CG_ABS_TOL,
	      const double rel_err = CG_REL_TOL,
	      const int max_its = CG_MAX_ITS);

void ell_mvp(const ell_matrix *m, const double *x, double *y);
int ell_solve_cgpd(const ell_matrix *m, const double *b, double *x,
		   double *err_);
void ell_add_2D(ell_matrix *m, int ex, int ey, const double *Ae);
void ell_add_3D(ell_matrix *m, int ex, int ey, int ez, const double *Ae);
void ell_set_zero_mat(ell_matrix *m);
void ell_set_bc_2D(ell_matrix *m);
void ell_set_bc_3D(ell_matrix *m);
void ell_free(ell_matrix *m);

double get_norm(const double *vector, const int n);
double get_dot(const double *v1, const double *v2, const int n);
double ell_get_norm(const ell_matrix *m);

void print_ell(const ell_matrix *A);

// OpenACC compatibility

void ell_mvp_acc(const ell_matrix *m, const double *x, double *y);

int ell_solve_cgpd_acc(const ell_matrix *m, const double *b, double *x,
		       double *err_);

double get_norm_acc(const double *vector, const int n);

double get_dot_acc(const double *v1, const double *v2, const int n);

void print_ell_acc(const ell_matrix *A);

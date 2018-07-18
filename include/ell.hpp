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

#ifndef ELL_H_
#define ELL_H_

#define nod_index(i,j,k) ((k)*nx*ny + (j)*nx + (i))

typedef struct {
	int nrow;		// number of rows
	int ncol;		// number of columns
	int nnz;		// non zeros per row
	int *cols;
	double *vals;
} ell_matrix;

typedef struct {
	int max_its;
	int its;
	double min_tol;
	double err;
	double *k, *r, *z, *p, *q;
} ell_solver;

int ell_add_val(ell_matrix *m, int row, int col, double val);
int ell_add_vals(ell_matrix *m, int *ix, int nx, int *iy, int ny, double *vals);

void ell_mvp(const ell_matrix *m, const double *x, double *y);
int ell_get_val(ell_matrix *m, int row, int col, double *val);

int ell_solve_cg(ell_solver *solver, ell_matrix *m, double *b, double *x);
int ell_solve_jacobi(ell_solver *solver, ell_matrix *m, double *b, double *x);
int ell_solve_jacobi_2D(ell_solver *solver, ell_matrix *m, int nFields,
                        int nx, int ny, double *b, double *x);

int ell_set_val(ell_matrix *m, int row, int col, double val);
int ell_set_zero_row(ell_matrix *m, int row, double diag_val);
int ell_set_zero_col(ell_matrix *m, int col, double diag_val);
int ell_set_zero_mat(ell_matrix *m);

int ell_print(ell_matrix *m);
int ell_print_full(ell_matrix *m);

void ell_free(ell_matrix *m);

int ell_solve_cgpd_struct(ell_solver *solver, const ell_matrix * const m,
                          int nFields, int dim, int nn,
                          const double * const b, double *x);

void ell_add_struct2D(ell_matrix *m, int ex, int ey, double *Ae,
                      int nFields, int nx, int ny);

void ell_add_struct3D(ell_matrix *m, int ex, int ey, int ez, double *Ae,
                      int nFields, int nx, int ny, int nz);

void ell_set_bc_2D(ell_matrix *m, int nFields, int nx, int ny);
void ell_set_bc_3D(ell_matrix *m, int nFields, int nx, int ny, int nz);

void ell_init_2D(ell_matrix *m, int nFields, int nx, int ny);
void ell_init_3D(ell_matrix *m, int nFields, int nx, int ny, int nz);

#endif

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

#include <iostream>
#include <iomanip>		// print with format

#include <cmath>

#include "ell.hpp"

#define nod_index(i,j,k) ((k)*nx*ny + (j)*nx + (i))

using namespace std;

void ell_free(ell_matrix &m)
{
	if (m.cols != NULL)
		free(m.cols);
	if (m.vals != NULL)
		free(m.vals);
}

int ell_set_zero_mat(ell_matrix *m)
{
	for (int i = 0; i < m->nrow; i++)
		for (int j = 0; j < m->nnz; j++)
			m->vals[i * m->nnz + j] = 0.0;
	return 0;
}

void ell_mvp_2D(ell_matrix *m, double *x, double *y)
{
	//  y = m * x
	for (int i = 0; i < m->nrow; i++) {
		y[i] = 0;
		for (int j = 0; j < m->nnz; j++)
			y[i] += m->vals[(i * m->nnz) + j] * x[m->cols[(i * m->nnz) + j]];
	}
}

int ell_solve_jacobi_2D(ell_solver *solver, ell_matrix *m, int nFields, int nx, int ny, double *b, double *x)
{
	/* A = K - N
	 * K = diag(A)
	 * N_ij = -a_ij for i!=j  and =0 if i=j
	 * x_(i) = K^-1 * ( N * x_(i-1) + b )
	 */
	if (m == NULL || b == NULL || x == NULL)
		return 1;

	double *k = (double *)malloc(m->nrow * sizeof(double));	// K = diag(A)
	double *e_i = (double *)malloc(m->nrow * sizeof(double));

	int nn = nx * ny;
	for (int i = 0; i < nn; i++) {
		for (int d = 0; d < nFields; d++) {
			k[i * nFields + d] = 1 / m->vals[i * nFields * m->nnz + 4 * nFields + d * m->nnz + d];
		}
	}

	int its = 0;
	int max_its = solver->max_its;
	double err;
	double min_tol = solver->min_tol;

	while (its < max_its) {
		err = 0;
		int i = 0;
		while (i < m->nrow) {
			double aux = 0.0;	// sum_(j!=i) a_ij * x_j
			int j = 0;
			while (j < m->nnz) {
				if (m->cols[i * m->nnz + j] == -1)
					break;
				if (m->cols[i * m->nnz + j] != i)
					aux += m->vals[i * m->nnz + j] * x[m->cols[i * m->nnz + j]];
				j++;
			}
			x[i] = k[i] * (-1 * aux + b[i]);
			i++;
		}

		err = 0;
		ell_mvp_2D(m, x, e_i);
		for (int i = 0; i < m->nrow; i++) {
			e_i[i] -= b[i];
			err += e_i[i] * e_i[i];
		}
		err = sqrt(err);
		if (err < min_tol)
			break;
		its++;
	}
	solver->err = err;
	solver->its = its;
	return 0;
}

int ell_solve_cgpd_2D(ell_solver *solver, ell_matrix *m, int nFields, int nx, int ny, double *b, double *x)
{
	/* cg with jacobi preconditioner
	 * r_1 residue in actual iteration
	 * z_1 = K^-1 * r_0 actual auxiliar vector
	 * rho_0 rho_1 = r_0^t * z_1 previous and actual iner products <r_i, K^-1, r_i>
	 * p_1 actual search direction
	 * q_1 = A*p_1 auxiliar vector
	 * d_1 = rho_0 / (p_1^t * q_1) actual step
	 * x_1 = x_0 - d_1 * p_1
	 * r_1 = r_0 - d_1 * q_1
	 */
	if (m == NULL || b == NULL || x == NULL)
		return 1;

	int its = 0;
	double *k = (double *) malloc(m->nrow * sizeof(double));	// K = diag(A)
	double *r = (double *) malloc(m->nrow * sizeof(double));
	double *z = (double *) malloc(m->nrow * sizeof(double));
	double *p = (double *) malloc(m->nrow * sizeof(double));
	double *q = (double *) malloc(m->nrow * sizeof(double));
	double rho_0, rho_1, d;
	double err;

	int nn = nx * ny;

	for (int i = 0; i < nn; i++) {
		for (int d = 0; d < nFields; d++) {
			k[i * nFields + d] = 1 / m->vals[i * nFields * m->nnz + 4 * nFields + d * m->nnz + d];
		}
	}

	ell_mvp_2D(m, x, r);

	for (int i = 0; i < m->nrow; i++)
		r[i] -= b[i];

	do {

		err = 0;
		for (int i = 0; i < m->nrow; i++)
			err += r[i] * r[i];
		//    cout << "cg_err = " << err << endl;
		err = sqrt(err);
		if (err < solver->min_tol)
			break;

		for (int i = 0; i < m->nrow; i++)
			z[i] = k[i] * r[i];

		rho_1 = 0.0;
		for (int i = 0; i < m->nrow; i++)
			rho_1 += r[i] * z[i];

		if (its == 0) {
			for (int i = 0; i < m->nrow; i++)
				p[i] = z[i];
		} else {
			double beta = rho_1 / rho_0;
			for (int i = 0; i < m->nrow; i++)
				p[i] = z[i] + beta * p[i];
		}

		ell_mvp_2D(m, p, q);
		double aux = 0;
		for (int i = 0; i < m->nrow; i++)
			aux += p[i] * q[i];
		d = rho_1 / aux;

		for (int i = 0; i < m->nrow; i++) {
			x[i] -= d * p[i];
			r[i] -= d * q[i];
		}

		rho_0 = rho_1;
		its++;

	} while (its < solver->max_its);

	solver->err = err;
	solver->its = its;

	free(k);
	free(r);
	free(z);
	free(p);
	free(q);

	return 0;
}

int ell_solve_cgpd_struct(ell_solver *solver, ell_matrix *m, int nFields, int dim, int nn, double *b, double *x)
{
	/* cg with jacobi preconditioner
	 * r_1 residue in actual iteration
	 * z_1 = K^-1 * r_0 actual auxiliar vector
	 * rho_0 rho_1 = r_0^t * z_1 previous and actual iner products <r_i, K^-1, r_i>
	 * p_1 actual search direction
	 * q_1 = A*p_1 auxiliar vector
	 * d_1 = rho_0 / (p_1^t * q_1) actual step
	 * x_1 = x_0 - d_1 * p_1
	 * r_1 = r_0 - d_1 * q_1
	 */
	if (m == NULL || b == NULL || x == NULL)
		return 1;

	int its = 0;
	double *k = (double *)malloc(m->nrow * sizeof(double));	// K = diag(A)
	double *r = (double *)malloc(m->nrow * sizeof(double));
	double *z = (double *)malloc(m->nrow * sizeof(double));
	double *p = (double *)malloc(m->nrow * sizeof(double));
	double *q = (double *)malloc(m->nrow * sizeof(double));
	double rho_0, rho_1, d;
	double err;

	if (dim == 2) {
		for (int i = 0; i < nn; i++) {
			for (int d = 0; d < nFields; d++) {
				k[i * nFields + d] = 1 / m->vals[i * nFields * m->nnz + 4 * nFields + d * m->nnz + d];
			}
		}
	} else if (dim == 3) {
		for (int i = 0; i < nn; i++) {
			for (int d = 0; d < nFields; d++) {
				k[i * nFields + d] = 1 / m->vals[i * nFields * m->nnz + 13 * nFields + d * m->nnz + d];
			}
		}
	}

	ell_mvp_2D(m, x, r);
	for (int i = 0; i < m->nrow; i++)
		r[i] -= b[i];

	do {

		err = 0;
		for (int i = 0; i < m->nrow; i++)
			err += r[i] * r[i];
		err = sqrt(err);
		if (err < solver->min_tol)
			break;

		for (int i = 0; i < m->nrow; i++)
			z[i] = k[i] * r[i];

		rho_1 = 0.0;
		for (int i = 0; i < m->nrow; i++)
			rho_1 += r[i] * z[i];

		if (its == 0) {
			for (int i = 0; i < m->nrow; i++)
				p[i] = z[i];
		} else {
			double beta = rho_1 / rho_0;
			for (int i = 0; i < m->nrow; i++)
				p[i] = z[i] + beta * p[i];
		}

		ell_mvp_2D(m, p, q);
		double aux = 0;
		for (int i = 0; i < m->nrow; i++)
			aux += p[i] * q[i];
		d = rho_1 / aux;

		for (int i = 0; i < m->nrow; i++) {
			x[i] -= d * p[i];
			r[i] -= d * q[i];
		}

		rho_0 = rho_1;
		its++;

	} while (its < solver->max_its);

	solver->err = err;
	solver->its = its;

	free(k);
	free(r);
	free(z);
	free(p);
	free(q);

	return 0;
}

int ell_print(ell_matrix *m)
{
	if (m == NULL)
		return 1;
	if (m->vals == NULL || m->cols == NULL)
		return 2;

	cout << "Cols = " << endl;
	for (int i = 0; i < m->nrow; i++) {
		cout << setw(7) << "row= " << i;
		for (int j = 0; j < m->nnz; j++) {
			cout << setw(7) << setprecision(4) << m->cols[i * m->nnz + j] << " ";
		}
		cout << endl;
	}

	cout << "Vals = " << endl;
	for (int i = 0; i < m->nrow; i++) {
		cout << setw(7) << "row= " << i;
		for (int j = 0; j < m->nnz; j++) {
			cout << setw(7) << setprecision(4) << m->vals[i * m->nnz + j] << " ";
		}
		cout << endl;
	}
	return 0;
}

void ell_add_struct(ell_matrix &m, int ex, int ey, double *Ae, int nFields, int nx, int ny)
{
	// assembly Ae in 2D structured grid representation
	// nFields : number of scalar components on each node

	const int npe = 4;
	int nnz = m.nnz;
	int cols_row_0[4] = { 4, 5, 8, 7 };
	int cols_row_1[4] = { 3, 4, 7, 6 };
	int cols_row_2[4] = { 0, 1, 4, 3 };
	int cols_row_3[4] = { 1, 2, 5, 4 };

	int n0 = ey * nx + ex;
	int n1 = ey * nx + ex + 1;
	int n2 = (ey + 1) * nx + ex + 1;
	int n3 = (ey + 1) * nx + ex;

	for (int i = 0; i < nFields; i++) {
		for (int n = 0; n < npe; n++) {
			for (int j = 0; j < nFields; j++) {
				m.vals[n0 * nFields * nnz + i * nnz + cols_row_0[n] * nFields + j] += Ae[0 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n1 * nFields * nnz + i * nnz + cols_row_1[n] * nFields + j] += Ae[1 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n2 * nFields * nnz + i * nnz + cols_row_2[n] * nFields + j] += Ae[2 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n3 * nFields * nnz + i * nnz + cols_row_3[n] * nFields + j] += Ae[3 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
			}
		}
	}

}

void ell_add_struct(ell_matrix &m, int ex, int ey, int ez, double *Ae, int nFields, int nx, int ny, int nz)
{
	// assembly Ae in 3D structured grid representation
	// nFields : number of scalar components on each node

	int npe = 8;

	int nnz = m.nnz;
	int cols_row_0[8] = { 13, 14, 17, 16, 22, 23, 26, 25 };
	int cols_row_1[8] = { 12, 13, 16, 15, 21, 22, 25, 24 };
	int cols_row_2[8] = { 9, 10, 13, 12, 18, 19, 22, 21 };
	int cols_row_3[8] = { 10, 11, 14, 13, 19, 20, 23, 22 };
	int cols_row_4[8] = { 4, 5, 8, 7, 13, 14, 17, 16 };
	int cols_row_5[8] = { 3, 4, 7, 6, 12, 13, 16, 15 };
	int cols_row_6[8] = { 0, 1, 4, 3, 9, 10, 13, 12 };
	int cols_row_7[8] = { 1, 2, 5, 4, 10, 11, 14, 13 };

	int n0 = ez * (nx * ny) + ey * nx + ex;
	int n1 = ez * (nx * ny) + ey * nx + ex + 1;
	int n2 = ez * (nx * ny) + (ey + 1) * nx + ex + 1;
	int n3 = ez * (nx * ny) + (ey + 1) * nx + ex;
	int n4 = n0 + nx * ny;
	int n5 = n1 + nx * ny;
	int n6 = n2 + nx * ny;
	int n7 = n3 + nx * ny;

	for (int i = 0; i < nFields; i++) {
		for (int n = 0; n < npe; n++) {
			for (int j = 0; j < nFields; j++) {
				m.vals[n0 * nFields * nnz + i * nnz + cols_row_0[n] * nFields + j] += Ae[0 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n1 * nFields * nnz + i * nnz + cols_row_1[n] * nFields + j] += Ae[1 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n2 * nFields * nnz + i * nnz + cols_row_2[n] * nFields + j] += Ae[2 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n3 * nFields * nnz + i * nnz + cols_row_3[n] * nFields + j] += Ae[3 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n4 * nFields * nnz + i * nnz + cols_row_4[n] * nFields + j] += Ae[4 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n5 * nFields * nnz + i * nnz + cols_row_5[n] * nFields + j] += Ae[5 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n6 * nFields * nnz + i * nnz + cols_row_6[n] * nFields + j] += Ae[6 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
				m.vals[n7 * nFields * nnz + i * nnz + cols_row_7[n] * nFields + j] += Ae[7 * (npe * nFields) * nFields + i * (npe * nFields) + n * nFields + j];
			}
		}
	}

}

void ell_set_bc_2D(ell_matrix &m, int nFields, int nx, int ny)
{
	// Sets 1 on the diagonal of the boundaries and does 0 on the columns corresponding to that values

	int nnz = m.nnz;

	// y=0
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < nx; n++) {
			for (int j = 0; j < nnz; j++) {
				m.vals[n * nFields * nnz + d * nnz + j] = 0;
			}
			m.vals[n * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}
	// y=ly
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < nx; n++) {
			for (int j = 0; j < nnz; j++) {
				m.vals[(n + (ny - 1) * nx) * nFields * nnz + d * nnz + j] = 0;
			}
			m.vals[(n + (ny - 1) * nx) * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}
	// x=0
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < ny - 2; n++) {
			for (int j = 0; j < nnz; j++) {
				m.vals[(n + 1) * nx * nFields * nnz + d * nnz + j] = 0;
			}
			m.vals[(n + 1) * nx * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}
	// x=lx
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < ny - 2; n++) {
			for (int j = 0; j < nnz; j++) {
				m.vals[((n + 2) * nx - 1) * nFields * nnz + d * nnz + j] = 0;
			}
			m.vals[((n + 2) * nx - 1) * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}

	// internal nodes next to the boundary

	// y = hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = nx + 2;
		for (int i = 0; i < nx - 4; i++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 1 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
			}
			n += 1;
		}
	}

	// y = ly - hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = (ny - 2) * nx + 2;
		for (int i = 0; i < nx - 4; i++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 7 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
			}
			n += 1;
		}
	}

	// x = hx
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = 2 * nx + 1;
		for (int i = 0; i < ny - 4; i++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 3 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
			}
			n += nx;
		}
	}

	// x = lx - hx
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = 4 * nx - 1;
		for (int i = 0; i < ny - 4; i++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 5 * nFields + d2] = 0;
				m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
			}
			n += nx;
		}
	}

	// x = hx , y = hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = nx + 1;
		for (int d2 = 0; d2 < nFields; d2++) {
			m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 1 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 3 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
		}
	}

	// x = lx - hx , y = hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = 2 * nx - 2;
		for (int d2 = 0; d2 < nFields; d2++) {
			m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 1 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 5 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
		}
	}

	// x = lx - hx , y = ly - hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = (ny - 1) * nx - 1;
		for (int d2 = 0; d2 < nFields; d2++) {
			m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 5 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 7 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
		}
	}

	// x = hx , y = ly - hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = (ny - 2) * nx + 1;
		for (int d2 = 0; d2 < nFields; d2++) {
			m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 3 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 7 * nFields + d2] = 0;
			m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
		}
	}

}

void ell_set_bc_3D(ell_matrix &m, int nFields, int nx, int ny, int nz)
{
	// Sets 1 on the diagonal of the boundaries and does 0 on the columns corresponding to that values
	int nnz = m.nnz;
	int n;

	// z=0
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				n = nod_index(i, j, 0);
				for (int col = 0; col < nnz; col++) {
					m.vals[n * nFields * nnz + d * nnz + col] = 0;
				}
				m.vals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// z=lz
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				n = nod_index(i, j, nz - 1);
				for (int col = 0; col < nnz; col++) {
					m.vals[n * nFields * nnz + d * nnz + col] = 0;
				}
				m.vals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// y=0
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(i, 0, k);
				for (int col = 0; col < nnz; col++) {
					m.vals[n * nFields * nnz + d * nnz + col] = 0;
				}
				m.vals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// y=ly
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(i, ny - 1, k);
				for (int col = 0; col < nnz; col++) {
					m.vals[n * nFields * nnz + d * nnz + col] = 0;
				}
				m.vals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// x=0
	for (int d = 0; d < nFields; d++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(0, j, k);
				for (int col = 0; col < nnz; col++) {
					m.vals[n * nFields * nnz + d * nnz + col] = 0;
				}
				m.vals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// x=lx
	for (int d = 0; d < nFields; d++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(nx - 1, j, k);
				for (int col = 0; col < nnz; col++) {
					m.vals[n * nFields * nnz + d * nnz + col] = 0;
				}
				m.vals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// z= 0 + dz
	if (nz > 2) {
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int i = 1; i < nx - 1; i++) {
				for (int j = 1; j < ny - 1; j++) {
					n = nod_index(i, j, 1);
					for (int d2 = 0; d2 < nFields; d2++) {
						m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 1 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 3 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 4 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 5 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 7 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
					}
				}
			}
		}
	}
	// z= lz - dz
	if (nz > 2) {
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int i = 1; i < nx - 1; i++) {
				for (int j = 1; j < ny - 1; j++) {
					n = nod_index(i, j, nz - 2);
					for (int d2 = 0; d2 < nFields; d2++) {
						m.vals[n * nFields * nnz + d1 * nnz + 18 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 19 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 20 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 21 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 22 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 23 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 24 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 25 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 26 * nFields + d2] = 0;
					}
				}
			}
		}
	}
	// y= 0 + dy
	if (ny > 2) {
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int i = 1; i < nx - 1; i++) {
				for (int k = 1; k < nz - 1; k++) {
					n = nod_index(i, 1, k);
					for (int d2 = 0; d2 < nFields; d2++) {
						m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 1 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 9 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 10 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 11 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 18 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 19 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 20 * nFields + d2] = 0;
					}
				}
			}
		}
	}
	// y= ly - dy
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int i = 1; i < nx - 1; i++) {
			for (int k = 1; k < nz - 1; k++) {
				n = nod_index(i, ny - 2, k);
				for (int d2 = 0; d2 < nFields; d2++) {
					m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 7 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 15 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 16 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 17 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 24 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 25 * nFields + d2] = 0;
					m.vals[n * nFields * nnz + d1 * nnz + 26 * nFields + d2] = 0;
				}
			}
		}
	}

	// x= 0 + dy
	if (nx > 2) {
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 1; k < nz - 1; k++) {
					n = nod_index(1, j, k);
					for (int d2 = 0; d2 < nFields; d2++) {
						m.vals[n * nFields * nnz + d1 * nnz + 0 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 3 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 9 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 12 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 15 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 18 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 21 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 24 * nFields + d2] = 0;
					}
				}
			}
		}
	}
	// x= lx - dx
	if (nx > 2) {
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 1; k < nz - 1; k++) {
					n = nod_index(nx - 2, j, k);
					for (int d2 = 0; d2 < nFields; d2++) {
						m.vals[n * nFields * nnz + d1 * nnz + 2 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 5 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 8 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 11 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 14 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 17 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 20 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 23 * nFields + d2] = 0;
						m.vals[n * nFields * nnz + d1 * nnz + 26 * nFields + d2] = 0;
					}
				}
			}
		}
	}

}

void ell_init_2D(ell_matrix &m, int nFields, int nx, int ny)
{
	int nn = nx * ny;
	m.nnz = 9 * nFields;
	m.nrow = nn * nFields;
	m.ncol = nn * nFields;
	m.cols = (int *)malloc(m.nnz * m.nrow * sizeof(int));
	m.vals = (double *)malloc(m.nnz * m.nrow * sizeof(double));

	int nnz = m.nnz;

	for (int i = 0; i < nn; i++) {
		// the coorners
		if (i == 0) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		} else if (i == nx - 1) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		} else if (i == nx * ny - 1) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		} else if (i == (ny - 1) * nx) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
		// y=0
		else if (i < nx) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		}
		// y=ly
		else if (i > (ny - 1) * nx) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
		// x=0
		else if ((i % nx) == 0) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		}
		// x=ly
		else if ((i + 1) % nx == 0) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
		// internal node
		else {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[i * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					m.cols[i * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					m.cols[i * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					m.cols[i * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					m.cols[i * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					m.cols[i * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		}
	}

}

void ell_init_3D(ell_matrix & m, int nFields, int nx, int ny, int nz)
{
	int nn = nx * ny * nz;
	m.nnz = 27 * nFields;
	m.nrow = nn * nFields;
	m.ncol = nn * nFields;
	m.cols = (int *)malloc(m.nnz * m.nrow * sizeof(int));
	m.vals = (double *)malloc(m.nnz * m.nrow * sizeof(double));

	int nnz = m.nnz;

	// x=0 y=0 z=0
	nn = nod_index(0, 0, 0);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
		}
	}

	// x=lx y=0 z=0
	nn = nod_index(nx - 1, 0, 0);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=lx y=ly z=0
	nn = nod_index(nx - 1, ny - 1, 0);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=ly z=0
	nn = nod_index(0, ny - 1, 0);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=0 z=lz
	nn = nod_index(0, 0, nz - 1);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=lx y=0 z=lz
	nn = nod_index(nx - 1, 0, nz - 1);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=lx y=ly z=lz
	nn = nod_index(nx - 1, ny - 1, nz - 1);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=ly z=lz
	nn = nod_index(0, ny - 1, nz - 1);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
			m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=0 (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(0, 0, k);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
			}
		}
	}

	// x=lx y=0 (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(nx - 1, 0, k);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=lx y=ly (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(nx - 1, ny - 1, k);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=0 y=ly (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(0, ny - 1, k);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// y=0 z=0 (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, 0, 0);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
			}
		}
	}

	// y=ly z=0 (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, ny - 1, 0);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// y=ly z=lz (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, ny - 1, nz - 1);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// y=0 z=lz (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, 0, nz - 1);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=0 z=0 (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(0, j, 0);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
			}
		}
	}

	// x=lx z=0 (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(nx - 1, j, 0);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=lx z=lz (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(nx - 1, j, nz - 1);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=0 z=lz (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(0, j, nz - 1);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
				m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// z=0
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			nn = nod_index(i, j, 0);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
				}
			}
		}
	}

	// z=lz
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			nn = nod_index(i, j, nz - 1);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
	}

	// y=0
	for (int i = 1; i < nx - 1; i++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(i, 0, k);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
				}
			}
		}
	}

	// y=ly
	for (int i = 1; i < nx - 1; i++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(i, ny - 1, k);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
	}

	// x=0
	for (int j = 1; j < ny - 1; j++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(0, j, k);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
				}
			}
		}
	}

	// x=lx
	for (int j = 1; j < ny - 1; j++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(nx - 1, j, k);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = 0;
					m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
					m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
	}

	// Internal nodes
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				nn = nod_index(i, j, k);
				for (int d1 = 0; d1 < nFields; d1++) {
					for (int d2 = 0; d2 < nFields; d2++) {
						m.cols[nn * nFields * nnz + 0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
						m.cols[nn * nFields * nnz + 2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
						m.cols[nn * nFields * nnz + 5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
						m.cols[nn * nFields * nnz + 8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
						m.cols[nn * nFields * nnz + 11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
						m.cols[nn * nFields * nnz + 14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
						m.cols[nn * nFields * nnz + 17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
						m.cols[nn * nFields * nnz + 20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
						m.cols[nn * nFields * nnz + 23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
						m.cols[nn * nFields * nnz + 25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
						m.cols[nn * nFields * nnz + 26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
					}
				}
			}
		}
	}
}

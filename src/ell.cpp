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

#include <iostream>
#include <iomanip>		// print with format

#include <cmath>

#include "ell.hpp"

#define nod_index(i,j,k) ((k)*nx*ny + (j)*nx + (i))

using namespace std;

void ell_free(ell_matrix *m)
{
	if (m->cols != NULL)
		free(m->cols);
	if (m->vals != NULL)
		free(m->vals);
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

int ell_solve_cgpd_struct(ell_solver *solver, ell_matrix *m, int nFields,
                          int dim, int nn, double *b, double *x)
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

void ell_add_struct2D(ell_matrix *m, int ex, int ey, double *Ae,
                    int nFields, int nx, int ny)
{
	// assembly Ae in 2D structured grid representation
	// nFields : number of scalar components on each node

	const int npe = 4;
	const int nnz = m->nnz;
	const int cols_row[4][4] = { { 4, 5, 8, 7 },
	                             { 3, 4, 7, 6 },
	                             { 0, 1, 4, 3 },
	                             { 1, 2, 5, 4 } };

	const int sn[4] = { ey * nx + ex,
	                   ey * nx + ex + 1,
	                   (ey + 1) * nx + ex + 1,
	                   (ey + 1) * nx + ex };

	const int nFieldsnnz = nFields * nnz;
	const int npenFields = npe * nFields;
	const int npenFields2 = npe * nFields * nFields;

	for (int i = 0; i < nFields; ++i)
		for (int n = 0; n < npe; ++n)
			for (int k = 0; k < 4; ++k)
				for (int j = 0; j < nFields; ++j)
					m->vals[sn[k] * nFieldsnnz + i * nnz + cols_row[k][n] * nFields + j]
						+= Ae[k * npenFields2 + i * npenFields + n * nFields + j];

}

void ell_add_struct3D(ell_matrix *m, int ex, int ey, int ez, double *Ae,
                    int nFields, int nx, int ny, int nz)
{
	// assembly Ae in 3D structured grid representation
	// nFields : number of scalar components on each node

	const int npe = 8;

	const int nnz = m->nnz;
	const int cols_row[8][8] = { { 13, 14, 17, 16, 22, 23, 26, 25 },
	                             { 12, 13, 16, 15, 21, 22, 25, 24 },
	                             { 9, 10, 13, 12, 18, 19, 22, 21 },
	                             { 10, 11, 14, 13, 19, 20, 23, 22 },
	                             { 4, 5, 8, 7, 13, 14, 17, 16 },
	                             { 3, 4, 7, 6, 12, 13, 16, 15 },
	                             { 0, 1, 4, 3, 9, 10, 13, 12 },
	                             { 1, 2, 5, 4, 10, 11, 14, 13 } };

	const int n0 = ez * (nx * ny) + ey * nx + ex;
	const int n1 = ez * (nx * ny) + ey * nx + ex + 1;
	const int n2 = ez * (nx * ny) + (ey + 1) * nx + ex + 1;
	const int n3 = ez * (nx * ny) + (ey + 1) * nx + ex;

	const int sn[8] = {	n0, n1, n2, n3,
		n0 + nx * ny,
		n1 + nx * ny,
		n2 + nx * ny,
		n3 + nx * ny };

	const int nFieldsnnz = nFields * nnz;
	const int npenFields = npe * nFields;
	const int npenFields2 = npe * nFields * nFields;

	for (int i = 0; i < nFields; ++i)
		for (int n = 0; n < npe; ++n)
			for (int k = 0; k < 8; ++k)
				for (int j = 0; j < nFields; ++j)
					m->vals[sn[k] * nFieldsnnz + i * nnz + cols_row[k][n] * nFields + j]
						+= Ae[k * npenFields2 + i * npenFields + n * nFields + j];
}

void ell_set_bc_2D(ell_matrix *m, int nFields, int nx, int ny)
{
	// Sets 1 on the diagonal of the boundaries and does 0 on the columns corresponding to that values

	const int nnz = m->nnz;
	double * const mvals = m->vals;

	// y=0
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < nx; n++) {
			for (int j = 0; j < nnz; j++) {
				mvals[n * nFields * nnz + d * nnz + j] = 0;
			}
			mvals[n * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}
	// y=ly
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < nx; n++) {
			for (int j = 0; j < nnz; j++) {
				mvals[(n + (ny - 1) * nx) * nFields * nnz + d * nnz + j] = 0;
			}
			mvals[(n + (ny - 1) * nx) * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}
	// x=0
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < ny - 2; n++) {
			for (int j = 0; j < nnz; j++) {
				mvals[(n + 1) * nx * nFields * nnz + d * nnz + j] = 0;
			}
			mvals[(n + 1) * nx * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}
	// x=lx
	for (int d = 0; d < nFields; d++) {
		for (int n = 0; n < ny - 2; n++) {
			for (int j = 0; j < nnz; j++) {
				mvals[((n + 2) * nx - 1) * nFields * nnz + d * nnz + j] = 0;
			}
			mvals[((n + 2) * nx - 1) * nFields * nnz + d * nnz + 4 * nFields + d] = 1;
		}
	}

	// internal nodes next to the boundary

	// y = hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = nx + 2;
		for (int i = 0; i < nx - 4; i++) {
			for (int l = 0; l < 3; ++l)
				for (int d2 = 0; d2 < nFields; d2++) 
					mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
			n += 1;
		}
	}

	// y = ly - hy
	for (int d1 = 0; d1 < nFields; d1++) {
		int n = (ny - 2) * nx + 2;
		for (int i = 0; i < nx - 4; i++) {
			for (int l = 6; l < 9; ++l)
				for (int d2 = 0; d2 < nFields; d2++)
					mvals[n * nFields * nnz + d1 * nnz + 6 * nFields + d2] = 0;
			n += 1;
		}
	}

	{ 	// x = hx
		const int tmp[] = { 0, 3, 6 };
		for (int d1 = 0; d1 < nFields; d1++) {
			int n = 2 * nx + 1;
			for (int i = 0; i < ny - 4; i++) {
				for (auto l : tmp)
					for (int d2 = 0; d2 < nFields; ++d2)
						mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
				n += nx;
			}
		}
	}

	{ // x = lx - hx
		const int tmp[] = { 2, 5, 8 };
		for (int d1 = 0; d1 < nFields; d1++) {
			int n = 4 * nx - 1;
			for (int i = 0; i < ny - 4; i++) {
				for (auto l : tmp)
					for (int d2 = 0; d2 < nFields; d2++)
						mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
				n += nx;
			}
		}
	}

	{ // x = hx , y = hy
		const int tmp[] = { 0, 1, 2, 3, 6 };
		for (int d1 = 0; d1 < nFields; d1++) {
			int n = nx + 1;
			for (auto l : tmp)
				for (int d2 = 0; d2 < nFields; d2++)
					mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;

		}
	}

	{// x = lx - hx , y = hy
		const int tmp[] = {0, 1, 2, 5, 8};
		for (int d1 = 0; d1 < nFields; d1++) {
			int n = 2 * nx - 2;
			for (auto l : tmp)
				for (int d2 = 0; d2 < nFields; d2++)
					mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
		}
	}

	{// x = lx - hx , y = ly - hy
		const int tmp[] = { 2, 5, 6, 7, 8 };
		for (int d1 = 0; d1 < nFields; d1++) {
			int n = (ny - 1) * nx - 1;
			for (auto l : tmp)
				for (int d2 = 0; d2 < nFields; d2++)
					mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
		}
	}

	{// x = hx , y = ly - hy
		const int tmp[] = { 0, 3, 6, 7, 8 };
		for (int d1 = 0; d1 < nFields; d1++) {
			int n = (ny - 2) * nx + 1;
			for (auto l : tmp)
				for (int d2 = 0; d2 < nFields; d2++)
					mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
		}
	}
}

void ell_set_bc_3D(ell_matrix *m, int nFields, int nx, int ny, int nz)
{
	// Sets 1 on the diagonal of the boundaries and does 0 on the columns corresponding to that values
	const int nnz = m->nnz;
	double * const mvals = m->vals;
	int n;

	// z=0
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				n = nod_index(i, j, 0);
				for (int col = 0; col < nnz; col++) {
					mvals[n * nFields * nnz + d * nnz + col] = 0;
				}
				mvals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// z=lz
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				n = nod_index(i, j, nz - 1);
				for (int col = 0; col < nnz; col++) {
					mvals[n * nFields * nnz + d * nnz + col] = 0;
				}
				mvals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// y=0
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(i, 0, k);
				for (int col = 0; col < nnz; col++) {
					mvals[n * nFields * nnz + d * nnz + col] = 0;
				}
				mvals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// y=ly
	for (int d = 0; d < nFields; d++) {
		for (int i = 0; i < nx; i++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(i, ny - 1, k);
				for (int col = 0; col < nnz; col++) {
					mvals[n * nFields * nnz + d * nnz + col] = 0;
				}
				mvals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// x=0
	for (int d = 0; d < nFields; d++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(0, j, k);
				for (int col = 0; col < nnz; col++) {
					mvals[n * nFields * nnz + d * nnz + col] = 0;
				}
				mvals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// x=lx
	for (int d = 0; d < nFields; d++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				n = nod_index(nx - 1, j, k);
				for (int col = 0; col < nnz; col++) {
					mvals[n * nFields * nnz + d * nnz + col] = 0;
				}
				mvals[n * nFields * nnz + d * nnz + 13 * nFields + d] = 1;
			}
		}
	}

	// z= 0 + dz
	if (nz > 2) {
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int i = 1; i < nx - 1; i++) {
				for (int j = 1; j < ny - 1; j++) {
					n = nod_index(i, j, 1);
					for (int l = 0; l < 9; ++l)
						for (int d2 = 0; d2 < nFields; d2++)
							mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
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
					for (int l = 18; l < 27; ++l)
						for (int d2 = 0; d2 < nFields; d2++)
							mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
				}
			}
		}
	}
	// y= 0 + dy
	if (ny > 2) {
		const int tmp[] = { 0, 1, 2, 9, 10, 11, 18, 19, 20 };
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int i = 1; i < nx - 1; i++) {
				for (int k = 1; k < nz - 1; k++) {
					n = nod_index(i, 1, k);
					for (auto l : tmp)
						for (int d2 = 0; d2 < nFields; d2++)
							mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
				}
			}
		}
	}
	// y= ly - dy
	{
		const int tmp[] = { 6, 7, 8, 15, 16, 17, 24, 25, 26 };
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int i = 1; i < nx - 1; i++) {
				for (int k = 1; k < nz - 1; k++) {
					n = nod_index(i, ny - 2, k);
					for (auto l : tmp)
						for (int d2 = 0; d2 < nFields; d2++)
							mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
				}
			}
		}
	}
	// x= 0 + dy
	if (nx > 2) {
		const int tmp[] = { 0, 3, 6, 9, 12, 15, 18, 21, 24 };
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 1; k < nz - 1; k++) {
					n = nod_index(1, j, k);
					for (auto l : tmp)
						for (int d2 = 0; d2 < nFields; ++d2)
							mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
				}
			}
		}
	}
	// x= lx - dx
	if (nx > 2) {
		const int tmp[] = { 2, 5, 8, 11, 14, 17, 20, 23, 26 };
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 1; k < nz - 1; k++) {
					n = nod_index(nx - 2, j, k);
					for (auto l : tmp )
						for (int d2 = 0; d2 < nFields; d2++)
							mvals[n * nFields * nnz + d1 * nnz + l * nFields + d2] = 0;
				}
			}
		}
	}

}

void ell_init_2D(ell_matrix *m, int nFields, int nx, int ny)
{
	const int nn = nx * ny;
	const int nnz = 9 * nFields;
	const int nrow = nn * nFields;

	m->nnz = nnz;
	m->nrow = nrow;
	m->ncol = nrow;
	m->cols = (int *) malloc(nnz * nrow * sizeof(int));
	m->vals = (double *) malloc(nnz * nrow * sizeof(double));


	for (int i = 0; i < nn; i++) {
		// the coorners
		int * const mcols = &(m->cols[i * nFields * nnz]);

		if (i == 0) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = 0;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = 0;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		} else if (i == nx - 1) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = 0;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = 0;
					mcols[6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					mcols[7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		} else if (i == nx * ny - 1) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					mcols[1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = 0;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = 0;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		} else if (i == (ny - 1) * nx) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					mcols[3 * nFields + nnz * d1 + d2] = 0;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = 0;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
		// y=0
		else if (i < nx) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = 0;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					mcols[7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		}
		// y=ly
		else if (i > (ny - 1) * nx) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					mcols[1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					mcols[3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = 0;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
		// x=0
		else if ((i % nx) == 0) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					mcols[3 * nFields + nnz * d1 + d2] = 0;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		}
		// x=ly
		else if ((i + 1) % nx == 0) {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					mcols[1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = 0;
					mcols[6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					mcols[7 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
		// internal node
		else {
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = (i - nx - 1) * nFields + d2;
					mcols[1 * nFields + nnz * d1 + d2] = (i - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = (i - nx + 1) * nFields + d2;
					mcols[3 * nFields + nnz * d1 + d2] = (i - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (i) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (i + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = (i + nx - 1) * nFields + d2;
					mcols[7 * nFields + nnz * d1 + d2] = (i + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = (i + nx + 1) * nFields + d2;
				}
			}
		}
	}

}

void ell_init_3D(ell_matrix *m, int nFields, int nx, int ny, int nz)
{
	int nn = nx * ny * nz;
	const int nnz = 27 * nFields;
	const int nrow = nn * nFields;
	int * mcols = NULL;

	m->nnz = nnz;
	m->nrow = nrow;
	m->ncol = nrow;
	m->cols = (int *) malloc(nnz * nrow * sizeof(int));
	m->vals = (double *) malloc(nnz * nrow * sizeof(double));

	// x=0 y=0 z=0
	nn = nod_index(0, 0, 0);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = 0;
			mcols[1 * nFields + nnz * d1 + d2] = 0;
			mcols[2 * nFields + nnz * d1 + d2] = 0;
			mcols[3 * nFields + nnz * d1 + d2] = 0;
			mcols[4 * nFields + nnz * d1 + d2] = 0;
			mcols[5 * nFields + nnz * d1 + d2] = 0;
			mcols[6 * nFields + nnz * d1 + d2] = 0;
			mcols[7 * nFields + nnz * d1 + d2] = 0;
			mcols[8 * nFields + nnz * d1 + d2] = 0;
			mcols[9 * nFields + nnz * d1 + d2] = 0;
			mcols[10 * nFields + nnz * d1 + d2] = 0;
			mcols[11 * nFields + nnz * d1 + d2] = 0;
			mcols[12 * nFields + nnz * d1 + d2] = 0;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			mcols[15 * nFields + nnz * d1 + d2] = 0;
			mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
			mcols[18 * nFields + nnz * d1 + d2] = 0;
			mcols[19 * nFields + nnz * d1 + d2] = 0;
			mcols[20 * nFields + nnz * d1 + d2] = 0;
			mcols[21 * nFields + nnz * d1 + d2] = 0;
			mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
			mcols[24 * nFields + nnz * d1 + d2] = 0;
			mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
			mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
		}
	}

	// x=lx y=0 z=0
	nn = nod_index(nx - 1, 0, 0);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = 0;
			mcols[1 * nFields + nnz * d1 + d2] = 0;
			mcols[2 * nFields + nnz * d1 + d2] = 0;
			mcols[3 * nFields + nnz * d1 + d2] = 0;
			mcols[4 * nFields + nnz * d1 + d2] = 0;
			mcols[5 * nFields + nnz * d1 + d2] = 0;
			mcols[6 * nFields + nnz * d1 + d2] = 0;
			mcols[7 * nFields + nnz * d1 + d2] = 0;
			mcols[8 * nFields + nnz * d1 + d2] = 0;
			mcols[9 * nFields + nnz * d1 + d2] = 0;
			mcols[10 * nFields + nnz * d1 + d2] = 0;
			mcols[11 * nFields + nnz * d1 + d2] = 0;
			mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = 0;
			mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
			mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			mcols[17 * nFields + nnz * d1 + d2] = 0;
			mcols[18 * nFields + nnz * d1 + d2] = 0;
			mcols[19 * nFields + nnz * d1 + d2] = 0;
			mcols[20 * nFields + nnz * d1 + d2] = 0;
			mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
			mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			mcols[23 * nFields + nnz * d1 + d2] = 0;
			mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
			mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
			mcols[26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=lx y=ly z=0
	nn = nod_index(nx - 1, ny - 1, 0);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = 0;
			mcols[1 * nFields + nnz * d1 + d2] = 0;
			mcols[2 * nFields + nnz * d1 + d2] = 0;
			mcols[3 * nFields + nnz * d1 + d2] = 0;
			mcols[4 * nFields + nnz * d1 + d2] = 0;
			mcols[5 * nFields + nnz * d1 + d2] = 0;
			mcols[6 * nFields + nnz * d1 + d2] = 0;
			mcols[7 * nFields + nnz * d1 + d2] = 0;
			mcols[8 * nFields + nnz * d1 + d2] = 0;
			mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
			mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			mcols[11 * nFields + nnz * d1 + d2] = 0;
			mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = 0;
			mcols[15 * nFields + nnz * d1 + d2] = 0;
			mcols[16 * nFields + nnz * d1 + d2] = 0;
			mcols[17 * nFields + nnz * d1 + d2] = 0;
			mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
			mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
			mcols[20 * nFields + nnz * d1 + d2] = 0;
			mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
			mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			mcols[23 * nFields + nnz * d1 + d2] = 0;
			mcols[24 * nFields + nnz * d1 + d2] = 0;
			mcols[25 * nFields + nnz * d1 + d2] = 0;
			mcols[26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=ly z=0
	nn = nod_index(0, ny - 1, 0);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = 0;
			mcols[1 * nFields + nnz * d1 + d2] = 0;
			mcols[2 * nFields + nnz * d1 + d2] = 0;
			mcols[3 * nFields + nnz * d1 + d2] = 0;
			mcols[4 * nFields + nnz * d1 + d2] = 0;
			mcols[5 * nFields + nnz * d1 + d2] = 0;
			mcols[6 * nFields + nnz * d1 + d2] = 0;
			mcols[7 * nFields + nnz * d1 + d2] = 0;
			mcols[8 * nFields + nnz * d1 + d2] = 0;
			mcols[9 * nFields + nnz * d1 + d2] = 0;
			mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
			mcols[12 * nFields + nnz * d1 + d2] = 0;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			mcols[15 * nFields + nnz * d1 + d2] = 0;
			mcols[16 * nFields + nnz * d1 + d2] = 0;
			mcols[17 * nFields + nnz * d1 + d2] = 0;
			mcols[18 * nFields + nnz * d1 + d2] = 0;
			mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
			mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
			mcols[21 * nFields + nnz * d1 + d2] = 0;
			mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
			mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
			mcols[24 * nFields + nnz * d1 + d2] = 0;
			mcols[25 * nFields + nnz * d1 + d2] = 0;
			mcols[26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=0 z=lz
	nn = nod_index(0, 0, nz - 1);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = 0;
			mcols[1 * nFields + nnz * d1 + d2] = 0;
			mcols[2 * nFields + nnz * d1 + d2] = 0;
			mcols[3 * nFields + nnz * d1 + d2] = 0;
			mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
			mcols[6 * nFields + nnz * d1 + d2] = 0;
			mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
			mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
			mcols[9 * nFields + nnz * d1 + d2] = 0;
			mcols[10 * nFields + nnz * d1 + d2] = 0;
			mcols[11 * nFields + nnz * d1 + d2] = 0;
			mcols[12 * nFields + nnz * d1 + d2] = 0;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			mcols[15 * nFields + nnz * d1 + d2] = 0;
			mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
			mcols[18 * nFields + nnz * d1 + d2] = 0;
			mcols[19 * nFields + nnz * d1 + d2] = 0;
			mcols[20 * nFields + nnz * d1 + d2] = 0;
			mcols[21 * nFields + nnz * d1 + d2] = 0;
			mcols[22 * nFields + nnz * d1 + d2] = 0;
			mcols[23 * nFields + nnz * d1 + d2] = 0;
			mcols[24 * nFields + nnz * d1 + d2] = 0;
			mcols[25 * nFields + nnz * d1 + d2] = 0;
			mcols[26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=lx y=0 z=lz
	nn = nod_index(nx - 1, 0, nz - 1);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = 0;
			mcols[1 * nFields + nnz * d1 + d2] = 0;
			mcols[2 * nFields + nnz * d1 + d2] = 0;
			mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
			mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			mcols[5 * nFields + nnz * d1 + d2] = 0;
			mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
			mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
			mcols[8 * nFields + nnz * d1 + d2] = 0;
			mcols[9 * nFields + nnz * d1 + d2] = 0;
			mcols[10 * nFields + nnz * d1 + d2] = 0;
			mcols[11 * nFields + nnz * d1 + d2] = 0;
			mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = 0;
			mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
			mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
			mcols[17 * nFields + nnz * d1 + d2] = 0;
			mcols[18 * nFields + nnz * d1 + d2] = 0;
			mcols[19 * nFields + nnz * d1 + d2] = 0;
			mcols[20 * nFields + nnz * d1 + d2] = 0;
			mcols[21 * nFields + nnz * d1 + d2] = 0;
			mcols[22 * nFields + nnz * d1 + d2] = 0;
			mcols[23 * nFields + nnz * d1 + d2] = 0;
			mcols[24 * nFields + nnz * d1 + d2] = 0;
			mcols[25 * nFields + nnz * d1 + d2] = 0;
			mcols[26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=lx y=ly z=lz
	nn = nod_index(nx - 1, ny - 1, nz - 1);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
			mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
			mcols[2 * nFields + nnz * d1 + d2] = 0;
			mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
			mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			mcols[5 * nFields + nnz * d1 + d2] = 0;
			mcols[6 * nFields + nnz * d1 + d2] = 0;
			mcols[7 * nFields + nnz * d1 + d2] = 0;
			mcols[8 * nFields + nnz * d1 + d2] = 0;
			mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
			mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			mcols[11 * nFields + nnz * d1 + d2] = 0;
			mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = 0;
			mcols[15 * nFields + nnz * d1 + d2] = 0;
			mcols[16 * nFields + nnz * d1 + d2] = 0;
			mcols[17 * nFields + nnz * d1 + d2] = 0;
			mcols[18 * nFields + nnz * d1 + d2] = 0;
			mcols[19 * nFields + nnz * d1 + d2] = 0;
			mcols[20 * nFields + nnz * d1 + d2] = 0;
			mcols[21 * nFields + nnz * d1 + d2] = 0;
			mcols[22 * nFields + nnz * d1 + d2] = 0;
			mcols[23 * nFields + nnz * d1 + d2] = 0;
			mcols[24 * nFields + nnz * d1 + d2] = 0;
			mcols[25 * nFields + nnz * d1 + d2] = 0;
			mcols[26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=ly z=lz
	nn = nod_index(0, ny - 1, nz - 1);
	mcols = &(m->cols[nn * nFields * nnz]);
	for (int d1 = 0; d1 < nFields; d1++) {
		for (int d2 = 0; d2 < nFields; d2++) {
			mcols[0 * nFields + nnz * d1 + d2] = 0;
			mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
			mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
			mcols[3 * nFields + nnz * d1 + d2] = 0;
			mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
			mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
			mcols[6 * nFields + nnz * d1 + d2] = 0;
			mcols[7 * nFields + nnz * d1 + d2] = 0;
			mcols[8 * nFields + nnz * d1 + d2] = 0;
			mcols[9 * nFields + nnz * d1 + d2] = 0;
			mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
			mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
			mcols[12 * nFields + nnz * d1 + d2] = 0;
			mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
			mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
			mcols[15 * nFields + nnz * d1 + d2] = 0;
			mcols[16 * nFields + nnz * d1 + d2] = 0;
			mcols[17 * nFields + nnz * d1 + d2] = 0;
			mcols[18 * nFields + nnz * d1 + d2] = 0;
			mcols[19 * nFields + nnz * d1 + d2] = 0;
			mcols[20 * nFields + nnz * d1 + d2] = 0;
			mcols[21 * nFields + nnz * d1 + d2] = 0;
			mcols[22 * nFields + nnz * d1 + d2] = 0;
			mcols[23 * nFields + nnz * d1 + d2] = 0;
			mcols[24 * nFields + nnz * d1 + d2] = 0;
			mcols[25 * nFields + nnz * d1 + d2] = 0;
			mcols[26 * nFields + nnz * d1 + d2] = 0;
		}
	}

	// x=0 y=0 (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(0, 0, k);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = 0;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = 0;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
				mcols[9 * nFields + nnz * d1 + d2] = 0;
				mcols[10 * nFields + nnz * d1 + d2] = 0;
				mcols[11 * nFields + nnz * d1 + d2] = 0;
				mcols[12 * nFields + nnz * d1 + d2] = 0;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = 0;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = 0;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = 0;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
			}
		}
	}

	// x=lx y=0 (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(nx - 1, 0, k);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = 0;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = 0;
				mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
				mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = 0;
				mcols[10 * nFields + nnz * d1 + d2] = 0;
				mcols[11 * nFields + nnz * d1 + d2] = 0;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = 0;
				mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = 0;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = 0;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = 0;
				mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
				mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=lx y=ly (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(nx - 1, ny - 1, k);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
				mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = 0;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = 0;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = 0;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = 0;
				mcols[15 * nFields + nnz * d1 + d2] = 0;
				mcols[16 * nFields + nnz * d1 + d2] = 0;
				mcols[17 * nFields + nnz * d1 + d2] = 0;
				mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
				mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = 0;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = 0;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=0 y=ly (linea)
	for (int k = 1; k < nz - 1; k++) {
		nn = nod_index(0, ny - 1, k);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
				mcols[3 * nFields + nnz * d1 + d2] = 0;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = 0;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = 0;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				mcols[12 * nFields + nnz * d1 + d2] = 0;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = 0;
				mcols[16 * nFields + nnz * d1 + d2] = 0;
				mcols[17 * nFields + nnz * d1 + d2] = 0;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
				mcols[21 * nFields + nnz * d1 + d2] = 0;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = 0;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// y=0 z=0 (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, 0, 0);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = 0;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = 0;
				mcols[4 * nFields + nnz * d1 + d2] = 0;
				mcols[5 * nFields + nnz * d1 + d2] = 0;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = 0;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = 0;
				mcols[10 * nFields + nnz * d1 + d2] = 0;
				mcols[11 * nFields + nnz * d1 + d2] = 0;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = 0;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
				mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
			}
		}
	}

	// y=ly z=0 (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, ny - 1, 0);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = 0;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = 0;
				mcols[4 * nFields + nnz * d1 + d2] = 0;
				mcols[5 * nFields + nnz * d1 + d2] = 0;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = 0;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = 0;
				mcols[16 * nFields + nnz * d1 + d2] = 0;
				mcols[17 * nFields + nnz * d1 + d2] = 0;
				mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
				mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
				mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = 0;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// y=ly z=lz (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, ny - 1, nz - 1);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
				mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
				mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = 0;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = 0;
				mcols[16 * nFields + nnz * d1 + d2] = 0;
				mcols[17 * nFields + nnz * d1 + d2] = 0;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = 0;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = 0;
				mcols[22 * nFields + nnz * d1 + d2] = 0;
				mcols[23 * nFields + nnz * d1 + d2] = 0;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = 0;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// y=0 z=lz (linea)
	for (int i = 1; i < nx - 1; i++) {
		nn = nod_index(i, 0, nz - 1);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = 0;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
				mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
				mcols[9 * nFields + nnz * d1 + d2] = 0;
				mcols[10 * nFields + nnz * d1 + d2] = 0;
				mcols[11 * nFields + nnz * d1 + d2] = 0;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = 0;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = 0;
				mcols[22 * nFields + nnz * d1 + d2] = 0;
				mcols[23 * nFields + nnz * d1 + d2] = 0;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = 0;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=0 z=0 (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(0, j, 0);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = 0;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = 0;
				mcols[4 * nFields + nnz * d1 + d2] = 0;
				mcols[5 * nFields + nnz * d1 + d2] = 0;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = 0;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = 0;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				mcols[12 * nFields + nnz * d1 + d2] = 0;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = 0;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
				mcols[21 * nFields + nnz * d1 + d2] = 0;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
			}
		}
	}

	// x=lx z=0 (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(nx - 1, j, 0);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = 0;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = 0;
				mcols[4 * nFields + nnz * d1 + d2] = 0;
				mcols[5 * nFields + nnz * d1 + d2] = 0;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = 0;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = 0;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = 0;
				mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = 0;
				mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
				mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
				mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
				mcols[23 * nFields + nnz * d1 + d2] = 0;
				mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
				mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=lx z=lz (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(nx - 1, j, nz - 1);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
				mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				mcols[2 * nFields + nnz * d1 + d2] = 0;
				mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = 0;
				mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
				mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				mcols[8 * nFields + nnz * d1 + d2] = 0;
				mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = 0;
				mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = 0;
				mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = 0;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = 0;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = 0;
				mcols[22 * nFields + nnz * d1 + d2] = 0;
				mcols[23 * nFields + nnz * d1 + d2] = 0;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = 0;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// x=0 z=lz (linea)
	for (int j = 1; j < ny - 1; j++) {
		nn = nod_index(0, j, nz - 1);
		mcols = &(m->cols[nn * nFields * nnz]);
		for (int d1 = 0; d1 < nFields; d1++) {
			for (int d2 = 0; d2 < nFields; d2++) {
				mcols[0 * nFields + nnz * d1 + d2] = 0;
				mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
				mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
				mcols[3 * nFields + nnz * d1 + d2] = 0;
				mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
				mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
				mcols[6 * nFields + nnz * d1 + d2] = 0;
				mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
				mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
				mcols[9 * nFields + nnz * d1 + d2] = 0;
				mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
				mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
				mcols[12 * nFields + nnz * d1 + d2] = 0;
				mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
				mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
				mcols[15 * nFields + nnz * d1 + d2] = 0;
				mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
				mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
				mcols[18 * nFields + nnz * d1 + d2] = 0;
				mcols[19 * nFields + nnz * d1 + d2] = 0;
				mcols[20 * nFields + nnz * d1 + d2] = 0;
				mcols[21 * nFields + nnz * d1 + d2] = 0;
				mcols[22 * nFields + nnz * d1 + d2] = 0;
				mcols[23 * nFields + nnz * d1 + d2] = 0;
				mcols[24 * nFields + nnz * d1 + d2] = 0;
				mcols[25 * nFields + nnz * d1 + d2] = 0;
				mcols[26 * nFields + nnz * d1 + d2] = 0;
			}
		}
	}

	// z=0
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			nn = nod_index(i, j, 0);
			mcols = &(m->cols[nn * nFields * nnz]);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = 0;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = 0;
					mcols[4 * nFields + nnz * d1 + d2] = 0;
					mcols[5 * nFields + nnz * d1 + d2] = 0;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = 0;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
					mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
					mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
					mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
					mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
				}
			}
		}
	}

	// z=lz
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			nn = nod_index(i, j, nz - 1);
			mcols = &(m->cols[nn * nFields * nnz]);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
					mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
					mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
					mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
					mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					mcols[18 * nFields + nnz * d1 + d2] = 0;
					mcols[19 * nFields + nnz * d1 + d2] = 0;
					mcols[20 * nFields + nnz * d1 + d2] = 0;
					mcols[21 * nFields + nnz * d1 + d2] = 0;
					mcols[22 * nFields + nnz * d1 + d2] = 0;
					mcols[23 * nFields + nnz * d1 + d2] = 0;
					mcols[24 * nFields + nnz * d1 + d2] = 0;
					mcols[25 * nFields + nnz * d1 + d2] = 0;
					mcols[26 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
	}

	// y=0
	for (int i = 1; i < nx - 1; i++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(i, 0, k);
			mcols = &(m->cols[nn * nFields * nnz]);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = 0;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
					mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
					mcols[9 * nFields + nnz * d1 + d2] = 0;
					mcols[10 * nFields + nnz * d1 + d2] = 0;
					mcols[11 * nFields + nnz * d1 + d2] = 0;
					mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					mcols[18 * nFields + nnz * d1 + d2] = 0;
					mcols[19 * nFields + nnz * d1 + d2] = 0;
					mcols[20 * nFields + nnz * d1 + d2] = 0;
					mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
					mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
				}
			}
		}
	}

	// y=ly
	for (int i = 1; i < nx - 1; i++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(i, ny - 1, k);
			mcols = &(m->cols[nn * nFields * nnz]);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
					mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
					mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = 0;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
					mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					mcols[15 * nFields + nnz * d1 + d2] = 0;
					mcols[16 * nFields + nnz * d1 + d2] = 0;
					mcols[17 * nFields + nnz * d1 + d2] = 0;
					mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
					mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
					mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					mcols[24 * nFields + nnz * d1 + d2] = 0;
					mcols[25 * nFields + nnz * d1 + d2] = 0;
					mcols[26 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
	}

	// x=0
	for (int j = 1; j < ny - 1; j++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(0, j, k);
			mcols = &(m->cols[nn * nFields * nnz]);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = 0;
					mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
					mcols[3 * nFields + nnz * d1 + d2] = 0;
					mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
					mcols[6 * nFields + nnz * d1 + d2] = 0;
					mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
					mcols[9 * nFields + nnz * d1 + d2] = 0;
					mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
					mcols[12 * nFields + nnz * d1 + d2] = 0;
					mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
					mcols[15 * nFields + nnz * d1 + d2] = 0;
					mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
					mcols[18 * nFields + nnz * d1 + d2] = 0;
					mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
					mcols[21 * nFields + nnz * d1 + d2] = 0;
					mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
					mcols[24 * nFields + nnz * d1 + d2] = 0;
					mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
				}
			}
		}
	}

	// x=lx
	for (int j = 1; j < ny - 1; j++) {
		for (int k = 1; k < nz - 1; k++) {
			nn = nod_index(nx - 1, j, k);
			mcols = &(m->cols[nn * nFields * nnz]);
			for (int d1 = 0; d1 < nFields; d1++) {
				for (int d2 = 0; d2 < nFields; d2++) {
					mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
					mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
					mcols[2 * nFields + nnz * d1 + d2] = 0;
					mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
					mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
					mcols[5 * nFields + nnz * d1 + d2] = 0;
					mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
					mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
					mcols[8 * nFields + nnz * d1 + d2] = 0;
					mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
					mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
					mcols[11 * nFields + nnz * d1 + d2] = 0;
					mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
					mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
					mcols[14 * nFields + nnz * d1 + d2] = 0;
					mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
					mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
					mcols[17 * nFields + nnz * d1 + d2] = 0;
					mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
					mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
					mcols[20 * nFields + nnz * d1 + d2] = 0;
					mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
					mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
					mcols[23 * nFields + nnz * d1 + d2] = 0;
					mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
					mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
					mcols[26 * nFields + nnz * d1 + d2] = 0;
				}
			}
		}
	}

	// Internal nodes
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				nn = nod_index(i, j, k);
				mcols = &(m->cols[nn * nFields * nnz]);
				for (int d1 = 0; d1 < nFields; d1++) {
					for (int d2 = 0; d2 < nFields; d2++) {
						mcols[0 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx - 1) * nFields + d2;
						mcols[1 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx) * nFields + d2;
						mcols[2 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - nx + 1) * nFields + d2;
						mcols[3 * nFields + nnz * d1 + d2] = (nn - (nx * ny) - 1) * nFields + d2;
						mcols[4 * nFields + nnz * d1 + d2] = (nn - (nx * ny)) * nFields + d2;
						mcols[5 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + 1) * nFields + d2;
						mcols[6 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx - 1) * nFields + d2;
						mcols[7 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx) * nFields + d2;
						mcols[8 * nFields + nnz * d1 + d2] = (nn - (nx * ny) + nx + 1) * nFields + d2;
						mcols[9 * nFields + nnz * d1 + d2] = (nn - nx - 1) * nFields + d2;
						mcols[10 * nFields + nnz * d1 + d2] = (nn - nx) * nFields + d2;
						mcols[11 * nFields + nnz * d1 + d2] = (nn - nx + 1) * nFields + d2;
						mcols[12 * nFields + nnz * d1 + d2] = (nn - 1) * nFields + d2;
						mcols[13 * nFields + nnz * d1 + d2] = (nn) * nFields + d2;
						mcols[14 * nFields + nnz * d1 + d2] = (nn + 1) * nFields + d2;
						mcols[15 * nFields + nnz * d1 + d2] = (nn + nx - 1) * nFields + d2;
						mcols[16 * nFields + nnz * d1 + d2] = (nn + nx) * nFields + d2;
						mcols[17 * nFields + nnz * d1 + d2] = (nn + nx + 1) * nFields + d2;
						mcols[18 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx - 1) * nFields + d2;
						mcols[19 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx) * nFields + d2;
						mcols[20 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - nx + 1) * nFields + d2;
						mcols[21 * nFields + nnz * d1 + d2] = (nn + (nx * ny) - 1) * nFields + d2;
						mcols[22 * nFields + nnz * d1 + d2] = (nn + (nx * ny)) * nFields + d2;
						mcols[23 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + 1) * nFields + d2;
						mcols[24 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx - 1) * nFields + d2;
						mcols[25 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx) * nFields + d2;
						mcols[26 * nFields + nnz * d1 + d2] = (nn + (nx * ny) + nx + 1) * nFields + d2;
					}
				}
			}
		}
	}
}

/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018
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


#include "ell.hpp"


int ell_solve_cgpd_cuda(const ell_matrix *m, const double *b, double *x, double *err)
{
	/* Conjugate Gradient Algorithm (CG) with Jacobi Preconditioner */

	if (!m || !b || !x)
		return 1;

	for (int i = 0; i < m->nn; i++) {
		for (int d = 0; d < m->nfield; d++)
			m->k[i * m->nfield + d] = 1 / m->vals[i * m->nfield * m->nnz
				+ m->shift * m->nfield + d * m->nnz + d];
	}

	for (int i = 0; i < m->nrow; ++i)
		x[i] = 0.0;

	ell_mvp(m, x, m->r);

	for (int i = 0; i < m->nrow; ++i)
		m->r[i] = b[i] - m->r[i];

	for (int i = 0; i < m->nrow; ++i)
		m->z[i] = m->k[i] * m->r[i];

	for (int i = 0; i < m->nrow; ++i)
		m->p[i] = m->z[i];

	double rz = get_dot(m->r, m->z, m->nrow);

	double pnorm_0 = sqrt(get_dot(m->z, m->z, m->nrow));
	double pnorm = pnorm_0;

	int its = 0;
	while (its < m->max_its) {

		if (pnorm < m->min_err || pnorm < pnorm_0 * m->rel_err)
			break;

		ell_mvp(m, m->p, m->Ap);
		double pAp = get_dot(m->p, m->Ap, m->nrow);

		const double alpha = rz / pAp;

		for (int i = 0; i < m->nrow; ++i)
			x[i] += alpha * m->p[i];

		for (int i = 0; i < m->nrow; ++i)
			m->r[i] -= alpha * m->Ap[i];

		for (int i = 0; i < m->nrow; ++i)
			m->z[i] = m->k[i] * m->r[i];

		pnorm = sqrt(get_dot(m->z, m->z, m->nrow));
		double rz_n = 0;
		rz_n = get_dot(m->r, m->z, m->nrow);

		const double beta = rz_n / rz;
		for (int i = 0; i < m->nrow; ++i)
			m->p[i] = m->z[i] + beta * m->p[i];

		rz = rz_n;
		its++;
	}

	*err = rz;

	return its;
}

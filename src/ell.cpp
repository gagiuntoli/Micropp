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

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "instrument.hpp"

using namespace std;

void ell_mvp(const ell_matrix *m, const double *x, double *y) {
  for (int i = 0; i < m->nrow; i++) {
    double tmp = 0;
    const int ix = i * m->nnz;
    for (int j = 0; j < m->nnz; j++) {
      tmp += m->vals[ix + j] * x[m->cols[ix + j]];
    }
    y[i] = tmp;
  }
}

double get_norm(const double *vector, const int n) {
  double norm = 0.0;
  for (int i = 0; i < n; ++i) norm += vector[i] * vector[i];
  return sqrt(norm);
}

double get_dot(const double *v1, const double *v2, const int n) {
  double prod = 0.0;
  for (int i = 0; i < n; ++i) prod += v1[i] * v2[i];
  return prod;
}

double ell_get_norm(const ell_matrix *m) {
  double norm = 0.0;
  for (int i = 0; i < m->nn; i++)
    for (int j = 0; j < m->nfield; j++)
      for (int d = 0; d < m->nnz; d++) norm += m->vals[i * m->nfield * m->nnz + j * m->nnz + d];
  return sqrt(norm);
}

int ell_solve_cgpd(const ell_matrix *m, const double *b, double *x, double *err) {
  INST_START;

  /* Conjugate Gradient Algorithm (CG) with Jacobi Preconditioner */

  if (!m || !b || !x) return 1;

  for (int i = 0; i < m->nn; i++) {
    for (int d = 0; d < m->nfield; d++)
      m->k[i * m->nfield + d] = 1 / m->vals[i * m->nfield * m->nnz + m->shift * m->nfield + d * m->nnz + d];
  }

  for (int i = 0; i < m->nrow; ++i) x[i] = 0.0;

  ell_mvp(m, x, m->r);

  for (int i = 0; i < m->nrow; ++i) m->r[i] = b[i] - m->r[i];

  for (int i = 0; i < m->nrow; ++i) m->z[i] = m->k[i] * m->r[i];

  for (int i = 0; i < m->nrow; ++i) m->p[i] = m->z[i];

  double rz = get_dot(m->r, m->z, m->nrow);

  double pnorm_0 = sqrt(get_dot(m->z, m->z, m->nrow));
  double pnorm = pnorm_0;

  int its = 0;
  while (its < m->max_its) {
    if (pnorm < m->min_err || pnorm < pnorm_0 * m->rel_err) break;

    ell_mvp(m, m->p, m->Ap);
    double pAp = get_dot(m->p, m->Ap, m->nrow);

    const double alpha = rz / pAp;

    for (int i = 0; i < m->nrow; ++i) x[i] += alpha * m->p[i];

    for (int i = 0; i < m->nrow; ++i) m->r[i] -= alpha * m->Ap[i];

    for (int i = 0; i < m->nrow; ++i) m->z[i] = m->k[i] * m->r[i];

    pnorm = sqrt(get_dot(m->z, m->z, m->nrow));

    double rz_n = get_dot(m->r, m->z, m->nrow);

    const double beta = rz_n / rz;
    for (int i = 0; i < m->nrow; ++i) m->p[i] = m->z[i] + beta * m->p[i];

    rz = rz_n;
    its++;
  }

  *err = rz;

  return its;
}

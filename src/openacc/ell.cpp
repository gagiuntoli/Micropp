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

#include "instrument.hpp"

using namespace std;

void ell_mvp_acc(const ell_matrix *m, const double *x, double *y) {
#pragma acc parallel loop gang worker present(m[ : 1], m->nrow, m->nnz, m->cols[ : m->nrow * m->nnz], \
                                              m->vals[ : m->nrow * m->nnz], x[ : m->nrow], y[ : m->nrow])
  for (int i = 0; i < m->nrow; i++) {
    double tmp = 0;
    const int ix = i * m->nnz;
#pragma acc loop vector
    for (int j = 0; j < m->nnz; j++) {
      tmp += m->vals[ix + j] * x[m->cols[ix + j]];
    }
    y[i] = tmp;
  }
}

double get_norm_acc(const double *vector, const int n) {
  double norm = 0.0;
#pragma acc parallel loop reduction(+ : norm) present(n, vector[ : n])
  for (int i = 0; i < n; ++i) norm += vector[i] * vector[i];
  return sqrt(norm);
}

double get_dot_acc(const double *v1, const double *v2, const int n) {
  double prod = 0.0;
#pragma acc parallel loop reduction(+ : prod) present(v1[ : n], v2[ : n])
  for (int i = 0; i < n; ++i) prod += v1[i] * v2[i];
  return prod;
}

int ell_solve_cgpd(const ell_matrix *m, const double *b, double *x, double *err) {
  INST_START;

  /* Conjugate Gradient Algorithm (CG) with Jacobi Preconditioner */

  if (!m || !b || !x) return 1;

#pragma acc enter data copyin(x[ : m->nrow], b[ : m->nrow])

#pragma acc enter data copyin(m[0 : 1])
#pragma acc enter data copyin(m->cols[ : m->nrow * m->nnz], m->vals[ : m->nrow * m->nnz])
#pragma acc enter data copyin(m->r[ : m->nrow], m->z[ : m->nrow], m->k[ : m->nrow], m->p[ : m->nrow], m->Ap[ : m->nrow])

#pragma acc parallel loop present(m[0 : 1], m->k[ : m->nrow], m->vals[ : m->nrow * m->nnz])
  for (int i = 0; i < m->nn; i++) {
    for (int d = 0; d < m->nfield; d++)
      m->k[i * m->nfield + d] = 1 / m->vals[i * m->nfield * m->nnz + m->shift * m->nfield + d * m->nnz + d];
  }

#pragma acc parallel loop present(m[0 : 1], m->nrow, x[ : m->nrow])
  for (int i = 0; i < m->nrow; ++i) x[i] = 0.0;

  ell_mvp_acc(m, x, m->r);

#pragma acc parallel loop present(m[0 : 1], b[ : m->nrow], m->r[m->nrow])
  for (int i = 0; i < m->nrow; ++i) m->r[i] = b[i] - m->r[i];

#pragma acc parallel loop present(m[0 : 1], m->z[ : m->nrow], m->k[ : m->nrow], m->r[m->nrow])
  for (int i = 0; i < m->nrow; ++i) m->z[i] = m->k[i] * m->r[i];

#pragma acc parallel loop present(m[0 : 1], m->p[ : m->nrow], m->z[m->nrow])
  for (int i = 0; i < m->nrow; ++i) m->p[i] = m->z[i];

  double rz = get_dot_acc(m->r, m->z, m->nrow);

  double pnorm_0 = sqrt(get_dot_acc(m->z, m->z, m->nrow));
  double pnorm = pnorm_0;

  int its = 0;
  while (its < m->max_its) {
    if (pnorm < m->min_err || pnorm < pnorm_0 * m->rel_err) break;

    ell_mvp_acc(m, m->p, m->Ap);
    double pAp = get_dot_acc(m->p, m->Ap, m->nrow);

    const double alpha = rz / pAp;

#pragma acc parallel loop present(m[0 : 1], x[ : m->nrow], m->p[m->nrow]) copyin(alpha)
    for (int i = 0; i < m->nrow; ++i) x[i] += alpha * m->p[i];

#pragma acc parallel loop present(m[0 : 1], m->r[ : m->nrow], m->Ap[m->nrow]) copyin(alpha)
    for (int i = 0; i < m->nrow; ++i) m->r[i] -= alpha * m->Ap[i];

#pragma acc parallel loop present(m[0 : 1], m->z[ : m->nrow], m->k[m->nrow], m->r[m->nrow])
    for (int i = 0; i < m->nrow; ++i) m->z[i] = m->k[i] * m->r[i];

    pnorm = sqrt(get_dot_acc(m->z, m->z, m->nrow));
    double rz_n = 0;
    rz_n = get_dot_acc(m->r, m->z, m->nrow);

    const double beta = rz_n / rz;
#pragma acc parallel loop present(m[0 : 1], m->z[ : m->nrow], m->p[m->nrow]) copyin(beta)
    for (int i = 0; i < m->nrow; ++i) m->p[i] = m->z[i] + beta * m->p[i];

    rz = rz_n;
    its++;
  }

#pragma acc exit data copyout(m->cols[ : m->nrow * m->nnz], m->vals[ : m->nrow * m->nnz])
#pragma acc exit data copyout(m->r[ : m->nrow], m->z[ : m->nrow], m->k[ : m->nrow], m->p[ : m->nrow], m->Ap[ : m->nrow])
#pragma acc exit data delete (m[0 : 1])

#pragma acc exit data copyout(x[ : m->nrow], b[ : m->nrow])

  *err = rz;

  return its;
}

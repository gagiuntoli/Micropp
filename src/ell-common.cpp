/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *						   Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "ell.hpp"
#include "instrument.hpp"

using namespace std;

void ell_init(ell_matrix *m, const int nfield, const int dim, const int ns[3], const double min_err,
              const double rel_err, const int max_its) {
  memcpy(m->n, ns, 3 * sizeof(int));
  assert(ns[0] >= 0 && ns[1] >= 0 && ns[2] >= 0);
  assert(dim >= 2 && dim <= 3);
  assert(nfield > 0);
  assert(max_its > 0);
  assert(min_err > 0);

  const int nx = ns[0];
  const int ny = ns[1];
  const int nz = ns[2];
  const int nn = (dim == 2) ? nx * ny : nx * ny * nz;
  const int nxny = nx * ny;
  const int num_nodes = (dim == 2) ? 9 : 27;
  const int nnz = num_nodes * nfield;
  const int nrow = nn * nfield;

  m->nn = nn;
  m->dim = dim;
  m->nfield = nfield;
  m->shift = (m->dim == 2) ? 4 : 13;
  m->nnz = nnz;
  m->nrow = nrow;
  m->ncol = nrow;
  m->cols = (int *)malloc(nnz * nrow * sizeof(int));
  m->vals = (double *)malloc(nnz * nrow * sizeof(double));

  m->max_its = max_its;
  m->min_err = min_err;
  m->rel_err = rel_err;
  m->k = (double *)malloc(nn * nfield * sizeof(double));
  m->r = (double *)malloc(nn * nfield * sizeof(double));
  m->z = (double *)malloc(nn * nfield * sizeof(double));
  m->p = (double *)malloc(nn * nfield * sizeof(double));
  m->Ap = (double *)malloc(nn * nfield * sizeof(double));

  if (dim == 2) {
    for (int fi = 0; fi < nfield; ++fi) {
      for (int xi = 0; xi < nx; ++xi) {
        for (int yi = 0; yi < ny; ++yi) {
          const int ni = nod_index2D(xi, yi);
          int *const cols_ptr = &(m->cols[ni * nfield * nnz + fi * nnz]);

          int ix[9] = {(yi == 0 || xi == 0) ? 0 : ni - nx - 1,
                       (yi == 0) ? 0 : ni - nx,
                       (yi == 0 || xi == nx - 1) ? 0 : ni - nx + 1,
                       (xi == 0) ? 0 : ni - 1,
                       ni,
                       (xi == nx - 1) ? 0 : ni + 1,
                       (xi == 0 || yi == ny - 1) ? 0 : ni + nx - 1,
                       (yi == ny - 1) ? 0 : ni + nx,
                       (xi == nx - 1 || yi == ny - 1) ? 0 : ni + nx + 1};

          for (int n = 0; n < num_nodes; ++n)
            for (int fj = 0; fj < nfield; ++fj) cols_ptr[n * nfield + fj] = ix[n] * nfield + fj;
        }
      }
    }

  } else if (dim == 3) {
    for (int fi = 0; fi < nfield; ++fi) {
      for (int xi = 0; xi < nx; ++xi) {
        for (int yi = 0; yi < ny; ++yi) {
          for (int zi = 0; zi < nz; ++zi) {
            const int ni = nod_index3D(xi, yi, zi);
            int *const cols_ptr = &(m->cols[ni * nfield * nnz + fi * nnz]);

            int ix[27] = {(zi == 0 || yi == 0 || xi == 0) ? 0 : ni - nxny - nx - 1,
                          (zi == 0 || yi == 0) ? 0 : ni - nxny - nx,
                          (zi == 0 || yi == 0 || xi == nx - 1) ? 0 : ni - nxny - nx + 1,
                          (zi == 0 || xi == 0) ? 0 : ni - nxny - 1,
                          (zi == 0) ? 0 : ni - nxny,
                          (zi == 0 || xi == nx - 1) ? 0 : ni - nxny + 1,
                          (zi == 0 || yi == ny - 1 || xi == 0) ? 0 : ni - nxny + nx - 1,
                          (zi == 0 || yi == ny - 1) ? 0 : ni - nxny + nx,
                          (zi == 0 || yi == ny - 1 || xi == nx - 1) ? 0 : ni - nxny + nx + 1,

                          (yi == 0 || xi == 0) ? 0 : ni - nx - 1,
                          (yi == 0) ? 0 : ni - nx,
                          (yi == 0 || xi == nx - 1) ? 0 : ni - nx + 1,
                          (xi == 0) ? 0 : ni - 1,
                          ni,
                          (xi == nx - 1) ? 0 : ni + 1,
                          (yi == ny - 1 || xi == 0) ? 0 : ni + nx - 1,
                          (yi == ny - 1) ? 0 : ni + nx,
                          (yi == ny - 1 || xi == nx - 1) ? 0 : ni + nx + 1,

                          (zi == nz - 1 || yi == 0 || xi == 0) ? 0 : ni + nxny - nx - 1,
                          (zi == nz - 1 || yi == 0) ? 0 : ni + nxny - nx,
                          (zi == nz - 1 || yi == 0 || xi == nx - 1) ? 0 : ni + nxny - nx + 1,
                          (zi == nz - 1 || xi == 0) ? 0 : ni + nxny - 1,
                          (zi == nz - 1) ? 0 : ni + nxny,
                          (zi == nz - 1 || xi == nx - 1) ? 0 : ni + nxny + 1,
                          (zi == nz - 1 || yi == ny - 1 || xi == 0) ? 0 : ni + nxny + nx - 1,
                          (zi == nz - 1 || yi == ny - 1) ? 0 : ni + nxny + nx,
                          (zi == nz - 1 || yi == ny - 1 || xi == nx - 1) ? 0 : ni + nxny + nx + 1};

            for (int n = 0; n < num_nodes; ++n)
              for (int fj = 0; fj < nfield; ++fj) cols_ptr[n * nfield + fj] = ix[n] * nfield + fj;
          }
        }
      }
    }
  }
}

void ell_add_2D(ell_matrix *m, int ex, int ey, const double *Ae) {
  // assembly Ae in 2D structured grid representation
  // nFields : number of scalar components on each node

  const int nx = m->n[0];
  const int nfield = m->nfield;
  const int npe = 4;
  const int nnz = m->nnz;
  const int cols_row[4][4] = {{4, 5, 8, 7}, {3, 4, 7, 6}, {0, 1, 4, 3}, {1, 2, 5, 4}};

  const int n0 = ey * nx + ex;
  const int ix_glo[4] = {n0, n0 + 1, n0 + nx + 1, n0 + nx};

  const int nnz_nfield = nfield * nnz;
  const int npe_nfield = npe * nfield;
  const int npe_nfield2 = npe * nfield * nfield;

  for (int fi = 0; fi < nfield; ++fi)
    for (int fj = 0; fj < nfield; ++fj)
      for (int i = 0; i < npe; ++i)
        for (int j = 0; j < npe; ++j)
          m->vals[ix_glo[i] * nnz_nfield + cols_row[i][j] * nfield + fi * nnz + fj] +=
              Ae[i * npe_nfield2 + fi * npe_nfield + j * nfield + fj];
}

void ell_add_3D(ell_matrix *m, int ex, int ey, int ez, const double *Ae) {
  // assembly Ae in 3D structured grid representation
  // nFields : number of scalar components on each node

  const int nx = m->n[0];
  const int ny = m->n[1];
  const int nfield = m->nfield;
  const int npe = 8;
  const int nnz = m->nnz;
  const int cols_row[8][8] = {{13, 14, 17, 16, 22, 23, 26, 25}, {12, 13, 16, 15, 21, 22, 25, 24},
                              {9, 10, 13, 12, 18, 19, 22, 21},  {10, 11, 14, 13, 19, 20, 23, 22},
                              {4, 5, 8, 7, 13, 14, 17, 16},     {3, 4, 7, 6, 12, 13, 16, 15},
                              {0, 1, 4, 3, 9, 10, 13, 12},      {1, 2, 5, 4, 10, 11, 14, 13}};

  const int nxny = nx * ny;
  const int n0 = ez * nxny + ey * nx + ex;
  const int n1 = n0 + 1;
  const int n2 = n0 + nx + 1;
  const int n3 = n0 + nx;

  const int ix_glo[8] = {n0, n1, n2, n3, n0 + nxny, n1 + nxny, n2 + nxny, n3 + nxny};

  const int nnz_nfield = nfield * nnz;
  const int npe_nfield = npe * nfield;
  const int npe_nfield2 = npe * nfield * nfield;

  for (int fi = 0; fi < nfield; ++fi)
    for (int fj = 0; fj < nfield; ++fj)
      for (int i = 0; i < npe; ++i)
        for (int j = 0; j < npe; ++j)
          m->vals[ix_glo[i] * nnz_nfield + cols_row[i][j] * nfield + fi * nnz + fj] +=
              Ae[i * npe_nfield2 + fi * npe_nfield + j * nfield + fj];
}

void ell_set_zero_mat(ell_matrix *m) { memset(m->vals, 0, m->nrow * m->nnz * sizeof(double)); }

void ell_set_bc_2D(ell_matrix *m) {
  // Sets 1s on the diagonal of the boundaries and 0s
  // on the columns corresponding to that values
  const int nx = m->n[0];
  const int ny = m->n[1];
  const int nfield = m->nfield;
  const int nnz = m->nnz;
  double *const mvals = m->vals;

  for (int d = 0; d < nfield; ++d) {
    for (int i = 0; i < nx; ++i) {
      const int n = nod_index2D(i, 0);  // y=0
      memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
      mvals[n * nfield * nnz + d * nnz + 4 * nfield + d] = 1;
    }

    for (int i = 0; i < nx; ++i) {
      const int n = nod_index2D(i, ny - 1);  // y=ly
      memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
      mvals[n * nfield * nnz + d * nnz + 4 * nfield + d] = 1;
    }

    for (int j = 1; j < ny - 1; ++j) {
      const int n = nod_index2D(0, j);  // x=0
      memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
      mvals[n * nfield * nnz + d * nnz + 4 * nfield + d] = 1;
    }

    for (int j = 1; j < ny - 1; ++j) {
      const int n = nod_index2D(nx - 1, j);  // x=lx
      memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
      mvals[n * nfield * nnz + d * nnz + 4 * nfield + d] = 1;
    }
  }
}

void ell_set_bc_3D(ell_matrix *m) {
  // Sets 1s on the diagonal of the boundaries and 0s
  // on the columns corresponding to that values
  const int nx = m->n[0];
  const int ny = m->n[1];
  const int nz = m->n[2];
  const int nfield = m->nfield;
  const int nnz = m->nnz;
  double *const mvals = m->vals;

  for (int d = 0; d < nfield; ++d) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        const int n = nod_index3D(i, j, 0);  // z=0
        memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
        mvals[n * nfield * nnz + d * nnz + 13 * nfield + d] = 1;
      }
    }

    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        const int n = nod_index3D(i, j, nz - 1);  // z=lz
        memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
        mvals[n * nfield * nnz + d * nnz + 13 * nfield + d] = 1;
      }
    }

    for (int i = 0; i < nx; ++i) {
      for (int k = 1; k < nz - 1; ++k) {
        const int n = nod_index3D(i, 0, k);  // y=0
        memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
        mvals[n * nfield * nnz + d * nnz + 13 * nfield + d] = 1;
      }
    }

    for (int i = 0; i < nx; ++i) {
      for (int k = 1; k < nz - 1; ++k) {
        const int n = nod_index3D(i, ny - 1, k);  // y=ly
        memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
        mvals[n * nfield * nnz + d * nnz + 13 * nfield + d] = 1;
      }
    }

    for (int j = 1; j < ny - 1; ++j) {
      for (int k = 1; k < nz - 1; ++k) {
        const int n = nod_index3D(0, j, k);  // x=0
        memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
        mvals[n * nfield * nnz + d * nnz + 13 * nfield + d] = 1;
      }
    }

    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        const int n = nod_index3D(nx - 1, j, k);  // x=lx
        memset(&mvals[n * nfield * nnz + d * nnz], 0, nnz * sizeof(double));
        mvals[n * nfield * nnz + d * nnz + 13 * nfield + d] = 1;
      }
    }
  }
}

void ell_free(ell_matrix *m) {
  if (m->cols != NULL) free(m->cols);
  if (m->vals != NULL) free(m->vals);
  if (m->k != NULL) free(m->k);
  if (m->r != NULL) free(m->r);
  if (m->z != NULL) free(m->z);
  if (m->p != NULL) free(m->p);
  if (m->Ap != NULL) free(m->Ap);
}

int ell_write(string filename, const ell_matrix *A) {
  ofstream file(filename, ios::out | ios::binary);
  if (!file) {
    cout << "Cannot open file:" << filename << endl;
    return 1;
  }
  file.write((char *)A, sizeof(ell_matrix));
  file.write((char *)A->vals, A->nrow * A->nnz * sizeof(double));
  file.write((char *)A->cols, A->nrow * A->nnz * sizeof(int));
  return 0;
}

int ell_read(string filename, ell_matrix *A) {
  ifstream file(filename, ios::in | ios::binary);
  if (!file) {
    cout << "Cannot open file:" << filename << endl;
    return 1;
  }
  file.read((char *)A, sizeof(ell_matrix));
  file.read((char *)A->vals, A->nrow * A->nnz * sizeof(double));
  file.read((char *)A->cols, A->nrow * A->nnz * sizeof(int));
  return 0;
}

void print_ell(const ell_matrix *A) {
  FILE *file;
  file = fopen("A.dat", "w");

  for (int i = 0; i < A->nrow; ++i)
    for (int j = 0; j < A->nnz; ++j)
      fprintf(file, "[%d][%d][%lf]\n", i, A->cols[i * A->nnz + j], A->vals[i * A->nnz + j]);

  fclose(file);
}

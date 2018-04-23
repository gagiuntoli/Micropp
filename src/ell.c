#include "ell.h"

int ell_init (ell_matrix * m, int nrow, int ncol, int nnz)
{
  if (m == NULL) return 1;
  m->nnz = nnz;
  m->nrow = nrow;
  m->ncol = ncol;
  m->cols = malloc((nrow*nnz) * sizeof(int));
  m->vals = malloc((nrow*nnz) * sizeof(double));
  if (m->vals == NULL || m->cols == NULL) return 2;
  for (int i = 0 ; i < (nrow*nnz) ; i++) m->cols[i] = -1;
  for (int i = 0 ; i < (nrow*nnz) ; i++) m->vals[i] = +0;
  return 0;
}

int ell_set_val (ell_matrix * m, int row, int col, double val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell error: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell error: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      m->cols[(row*m->nnz) + j] = col;
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == col) {
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    }
    j++;
  }
  if (j == m->nnz) {
    printf(RED "ell error: not enought space to store value in row %d and col %d\n" NRM, row, col);
    return 3;
  }
  return 4;
}

int ell_add_val (ell_matrix * m, int row, int col, double val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell error: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell error: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      m->cols[(row*m->nnz) + j] = col;
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == col) {
      m->vals[(row*m->nnz) + j] += val;
      return 0;
    }
    j++;
  }
  if (j == m->nnz) {
    printf(RED "ell error: not enought space to add value in row %d and col %d\n" NRM, row, col);
    return 3;
  }
  return 4;
}

int ell_add_vals (ell_matrix *m, int *ix, int nx, int *iy, int ny, double *vals)
{
  if (m == NULL || ix == NULL || iy == NULL || vals == NULL) return 1;
  for (int i = 0 ; i < nx ; i++) {
    for (int j = 0 ; j < ny ; j++) {
      ell_add_val (m, ix[i], iy[j], vals[i*nx + j]);
    }
  }
  return 0;
}

int ell_set_zero_row (ell_matrix *m, int row, double diag_val)
{
  if (m == NULL) return 1;
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == row) {
      m->vals[(row*m->nnz) + j] = diag_val;
    } else {
      m->vals[(row*m->nnz) + j] = 0.0;
    }
    j++;
  }
  return 0;
}

int ell_set_zero_col (ell_matrix *m, int col, double diag_val)
{
  if (m == NULL) return 1;
  for (int i = 0 ; i < m->nrow ; i++) {
    int j = 0;
    while (j < m->nnz) {
      if (m->cols[(i*m->nnz) + j] == -1) {
	break;
      } else if (m->cols[(i*m->nnz) + j] == col) {
	m->vals[(i*m->nnz) + j] = (i == col) ? diag_val : 0.0;
      }
      j++;
    }
  }
  return 0;
}

int ell_set_zero_mat (ell_matrix * m)
{
  for (int i = 0 ; i < m->nrow ; i++) ell_set_zero_row (m, i, 0.0);
  return 0;
}


int ell_mvp (ell_matrix * m, double *x, double *y)
{
  //  y = m * x
  if (m == NULL || x == NULL || y == NULL) return 1;

#pragma omp parallel for
  for (int i = 0 ; i < m->nrow ; i++) {
    y[i] = 0;
    int j = 0;
    while (j < m->nnz) {
      if (m->cols[(i*m->nnz) + j] == -1) break;
      y[i] += m->vals[(i*m->nnz) + j] * x[m->cols[(i*m->nnz) + j]];
      j++;
    }
  }
  return 0;
}

int ell_solve_jacobi (ell_solver *solver, ell_matrix * m, double *b, double *x)
{
  /* A = K - N  
   * K = diag(A)
   * N_ij = -a_ij for i!=j  and =0 if i=j
   * x_(i) = K^-1 * ( N * x_(i-1) + b )
   */
  if (m == NULL || b == NULL || x == NULL) return 1;
  
  double *k = malloc (m->nrow * sizeof(double)); // K = diag(A)
  double *e_i = malloc(m->nrow * sizeof(double));

  for (int i = 0 ; i < m->nrow ; i++) {
    ell_get_val (m, i, i, &k[i]);
    k[i] = 1 / k[i];
  }

  int its = 0;
  int max_its = solver->max_its;
  double err;
  double min_tol = solver->min_tol;

  while (its < max_its) {

    err = 0;
    int i = 0;
    while (i < m->nrow) {
      double aux = 0.0; // sum_(j!=i) a_ij * x_j
      int j = 0;
      while (j < m->nnz) {
        if (m->cols[i*m->nnz + j] == -1) break;
        if (m->cols[i*m->nnz + j] != i)
	  aux += m->vals[i*m->nnz + j] * x[m->cols[i*m->nnz + j]];
	j++;
      }
      x[i] = k[i] * (-1*aux + b[i]);
      i++;
    }

    err = 0;
    ell_mvp (m, x, e_i);
    for (int i = 0 ; i < m->nrow ; i++){
      e_i[i] -= b[i];
      err += e_i[i] * e_i[i];
    }
    err = sqrt(err); if (err < min_tol) break;
    its ++;
  }
  solver->err = err;
  solver->its = its;
  return 0;
}

int ell_solve_cg (ell_solver *solver, ell_matrix * m, double *b, double *x)
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
  if (m == NULL || b == NULL || x == NULL) return 1;
  
  int its = 0;
  double *k = malloc(m->nrow * sizeof(double)); // K = diag(A)
  double *r = malloc(m->nrow * sizeof(double));
  double *z = malloc(m->nrow * sizeof(double));
  double *p = malloc(m->nrow * sizeof(double));
  double *q = malloc(m->nrow * sizeof(double));
  double rho_0, rho_1, d;
  double err;

  for (int i = 0 ; i < m->nrow ; i++) {
    ell_get_val (m, i, i, &k[i]);
    k[i] = 1 / k[i];
  }

  ell_mvp (m, x, r);
  for (int i = 0 ; i < m->nrow ; i++)
    r[i] -= b[i];

  err = 0;
  for (int i = 0 ; i < m->nrow ; i++)
    err += r[i] * r[i];
  err = sqrt(err); if (err < solver->min_tol) return 0;

  while (its < solver->max_its) {

    for (int i = 0 ; i < m->nrow ; i++)
      z[i] = k[i] * r[i];

    rho_1 = 0.0;
    for (int i = 0 ; i < m->nrow ; i++)
      rho_1 += r[i] * z[i];

    if (its == 0) {
      for (int i = 0 ; i < m->nrow ; i++)
	p[i] = z[i];
    } else {
      double beta = rho_1 / rho_0;
      for (int i = 0 ; i < m->nrow ; i++)
	p[i] = z[i] + beta * p[i];
    }

    ell_mvp (m, p, q);
    double aux = 0;
    for (int i = 0 ; i < m->nrow ; i++)
      aux += p[i] * q[i];
    d = rho_1 / aux;

    for (int i = 0 ; i < m->nrow ; i++) {
      x[i] -= d * p[i];
      r[i] -= d * q[i];
    }

    rho_0 = rho_1;

    err = 0;
    for (int i = 0 ; i < m->nrow ; i++)
      err += r[i] * r[i];
    err = sqrt(err); if (err < solver->min_tol) break;

    its ++;
  }
  solver->err = err;
  solver->its = its;
  return 0;
}

int ell_get_val (ell_matrix * m, int row, int col, double *val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell_get_val: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell_get_val: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row * m->nnz) + j] == -1) {
      *val = 0.0;
      return 0;
    } else if (m->cols[(row * m->nnz) + j] == col) {
      *val = m->vals[(row * m->nnz) + j];
      return 0;
    }
    j++;
  }
  return 3;
}

int ell_print_full (ell_matrix * m)
{
  if (m == NULL) return 1;
  if (m->vals == NULL || m->cols == NULL) return 2;
  double val;
  for (int i = 0 ; i < m->nrow ; i++) {
    for (int j = 0 ; j < m->ncol ; j++) {
      printf("%lf%s",(ell_get_val(m, i, j, &val) == 0)?val:0.0, (j == m->ncol - 1) ? "\n" : " ");
    }
  }
  return 0;
}

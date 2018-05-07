#ifndef ELL_H_
#define ELL_H_

typedef struct {
  int nrow;     // number of rows
  int ncol;     // number of columns
  int nnz;      // non zeros per row
  int *cols;
  double *vals;
} ell_matrix;

typedef struct {
  int max_its;
  int its;
  double min_tol;
  double err;
} ell_solver;

int ell_set_zero_mat (ell_matrix * m);
int ell_print_full (ell_matrix * m);
int ell_print (ell_matrix * m);

int ell_init_2D (ell_matrix *m, int nrow, int ncol, int nnz);
void ell_mvp_2D (ell_matrix *m, double *x, double *y);
void ell_add_2D (ell_matrix &m, int e, double *Ae, int nFields, int nx, int ny);
void ell_init_2D (ell_matrix &m, int nFields, int nx, int ny);
void ell_set_bc_2D (ell_matrix &m, int nFields, int nx, int ny);
int ell_solve_cgpd_2D (ell_solver *solver, ell_matrix * m, int nFields, int nx, int ny, double *b, double *x);
int ell_solve_jacobi_2D (ell_solver *solver, ell_matrix * m, int nFields, int nx, int ny, double *b, double *x);

#endif

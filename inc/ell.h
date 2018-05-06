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

int ell_init (ell_matrix *m, int nrow, int ncol, int nnz);
int ell_set_val (ell_matrix *m, int row, int col, double val);
int ell_add_val (ell_matrix *m, int row, int col, double val);
int ell_add_vals (ell_matrix *m, int *ix, int nx, int *iy, int ny, double *vals);
int ell_mvp (ell_matrix *m, double *x, double *y);
int ell_get_val (ell_matrix *m, int row, int col, double *val);
int ell_solve_jacobi (ell_solver *solver, ell_matrix * m, double *b, double *x);
int ell_solve_cg (ell_solver *solver, ell_matrix * m, double *b, double *x);
int ell_set_zero_row (ell_matrix *m, int row, double diag_val);
int ell_set_zero_col (ell_matrix *m, int col, double diag_val);
int ell_set_zero_mat (ell_matrix * m);
int ell_print_full (ell_matrix * m);
int ell_print (ell_matrix * m);

void ell_add_2D (ell_matrix &m, int e, double *Ae, int nFields, int nx, int ny);
void ell_init_2D (ell_matrix &m, int nFields, int nx, int ny);

#endif

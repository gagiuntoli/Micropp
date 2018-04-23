#ifndef CSR_H_
#define CSR_H_

struct csr_matrix {
  unsigned int num_rows;
  unsigned int nnz;
  unsigned int *row_offsets;
  unsigned int *cols;
  double *coefs;
};

struct csr_vector {
  unsigned int n;
  double *coefs;
};

int csr_alloc_A (csr_matrix& A, int N);
int csr_free_A (csr_matrix& A);
void csr_alloc_v (csr_vector &v, unsigned int n);
void csr_free_v (csr_vector &v);
void csr_cg (const csr_matrix &A, const csr_vector& x, csr_vector &b);

#endif

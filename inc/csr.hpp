#include<cstdlib>
#include<cmath>

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

int csr_alloc_A (csr_matrix& A, int N)
{
  int num_rows=(N+1)*(N+1)*(N+1);
  int nnz=18*num_rows; // 9*2 = 18 if dim = 2 
  A.num_rows=num_rows;
  A.row_offsets=(unsigned int*)malloc((num_rows+1)*sizeof(unsigned int));
  A.cols=(unsigned int*)malloc(nnz*sizeof(unsigned int));
  A.coefs=(double*)malloc(nnz*sizeof(double));
  return 0;
}

int csr_free_A (csr_matrix& A)
{
  free(A.row_offsets);
  free(A.cols);
  free(A.coefs);
  return 0;
}

int get_elemental (int e, double (&Ae)[8][8])
{
  return 0;
}

int csr_assembly_A (csr_matrix& A)
{
  int dim = 3, npe = 4, nelem = 1; 
  int ix[8];
  double Ae[8][8];

  for (int e = 0; e<nelem; e++){
    get_elemental (e, Ae);
  }

  return 0;
}

void csr_alloc_v (csr_vector &v, unsigned int n)
{
  v.n=n;
  v.coefs=(double*)malloc(n*sizeof(double));
}

void csr_free_v (csr_vector &v)
{
  double *vcoefs=v.coefs;
  free(v.coefs);
}

void csr_mvp (const csr_matrix &A, const csr_vector& x, csr_vector &y)
{
  unsigned int num_rows=A.num_rows;
  unsigned int *restrict row_offsets=A.row_offsets;
  unsigned int *restrict cols=A.cols;
  double *restrict Acoefs=A.coefs;
  double *restrict xcoefs=x.coefs;
  double *restrict ycoefs=y.coefs;

  for(int i=0;i<num_rows;i++) {
    double sum=0;
    int row_start=row_offsets[i];
    int row_end=row_offsets[i+1];
    for(int j=row_start;j<row_end;j++) {
      unsigned int Acol=cols[j];
      double Acoef=Acoefs[j];
      double xcoef=xcoefs[Acol];
      sum+=Acoef*xcoef;
    }
    ycoefs[i]=sum;
  }
}

double dot(const csr_vector& x, const csr_vector &y)
{
  double sum=0;
  unsigned int n=x.n;
  double *restrict xcoefs=x.coefs;
  double *restrict ycoefs=y.coefs;

  for(int i=0;i<n;i++) {
    sum+=xcoefs[i]*ycoefs[i];
  }
  return sum;
}

void waxpby(double alpha, const csr_vector &x, double beta, const csr_vector &y, const csr_vector& w)
{
  // w = alpha * x + beta * y
  unsigned int n=x.n;
  double *restrict xcoefs=x.coefs;
  double *restrict ycoefs=y.coefs;
  double *restrict wcoefs=w.coefs;

  for(int i=0;i<n;i++) {
    wcoefs[i]=alpha*xcoefs[i]+beta*ycoefs[i];
  }
}

#define N 200
#define MAX_ITERS 100
#define TOL 1e-12

void csr_cg (const csr_matrix &A, const csr_vector& x, csr_vector &b)
{
  // solves A * x = b (finds x, initial guess also)
  double one=1.0, zero=0.0;
  double normr, rtrans, oldtrans, p_ap_dot , alpha, beta;
  int iter=0;

  waxpby(one, x, zero, x, p);
  matvec(A,p,Ap);
  waxpby(one, b, -one, Ap, r);
  
  rtrans=dot(r,r);
  normr=sqrt(rtrans);
  
  do {
    if(iter==0) {
      waxpby(one,r,zero,r,p);
    } else {
      oldtrans=rtrans;
      rtrans = dot(r,r);
      beta = rtrans/oldtrans;
      waxpby(one,r,beta,p,p);
    }
    
    normr=sqrt(rtrans);
  
    matvec(A,p,Ap);
    p_ap_dot = dot(Ap,p);

    alpha = rtrans/p_ap_dot;

    waxpby(one,x,alpha,p,x);
    waxpby(one,r,-alpha,Ap,r);

//    if(iter%10==0)
//      printf("Iteration: %d, Tolerance: %.4e\n", iter, normr);
    iter++;

  } while(iter<MAX_ITERS && normr>TOL);

  return 0;
}

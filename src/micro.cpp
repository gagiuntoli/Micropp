#include <vector>
#include <iostream>
#include "micro.h"

Problem::Problem (int dim, int size[3], int cg_its, double cg_tol)
{
  lx = 1.0; ly = 1.0; lz = 1.0;

  solver.max_its = cg_its;
  solver.min_tol = cg_tol;
  this->dim = dim;
  if (dim == 2) {
    nx = size[0];
    ny = size[1];
    nvoi = 3;
    nn = nx * ny;
    npe = 4;
  } else if (dim == 3) {
    nx = size[0];
    ny = size[1];
    nz = size[2];
    nvoi = 6;
    nn = nx * ny * nz;
    npe = 8;
  }

  dx = lx/(nx-1);
  dy = ly/(ny-1);
  if (dim == 2) {
    nelem = (nx-1) * (ny-1);
  } else if (dim == 3) {
    dz = lz/(nz-1);
    nelem = (nx-1) * (ny-1) * (nz-1);
  }

  b  = (double*)malloc(nn*dim*sizeof(double));
  du = (double*)malloc(nn*dim*sizeof(double));
  u  = (double*)malloc(nn*dim*sizeof(double));
  stress = (double*)malloc(nelem*nvoi*sizeof(double));
  strain = (double*)malloc(nelem*nvoi*sizeof(double));
  elem_type = (int*)malloc(nelem*sizeof(int));
//  int_vars  = (double*)malloc(nn*dim*sizeof(double));

  for (int i=0; i<nn*dim; i++)
    u[i] = 0.0;

  for (int e=0; e<nelem; e++)
    elem_type[e] = getElemType(e);

  flag_print_A = false;
  flag_print_b = false;
  flag_print_u = false;
  flag_print_du = false;
  flag_print_newton = false;
  flag_print_solver = false;

  NewRap_Its = 3;
  NewRap_Tol = 1.0e-5;

  if (dim == 2) {
    ell_init_2D (A, dim, nx, ny);
  } else if (dim == 3) {
    ell_init_3D (A, dim, nx, ny, nz);
  }

  return;
}

Problem::Problem (int dim, int size[3], int micro_type, double *micro_params, int *mat_types, double *params)
{
  lx = 1.0; ly = 1.0; lz = 1.0;

  this->micro_type = micro_type;

  int nParams;
  if (micro_type == 0) {
    // mat 1 = matrix
    // mat 2 = sphere
    numMaterials = 2;
    nParams = 4;
  }

  for (int i=0; i<nParams; i++) {
    this->micro_params[i] = micro_params[i];
  }
  lx = micro_params[0];
  ly = micro_params[1];
  lz = micro_params[2];

  for (int i=0; i<numMaterials; i++) {
    if (mat_types[i] == 0) {
      // lineal
      material_list[i].E  = params[i*MAX_MAT_PARAM + 0];
      material_list[i].nu = params[i*MAX_MAT_PARAM + 1];
      material_list[i].plasticity = false;
      material_list[i].damage     = false;
    } else if (mat_types[i] == 1) {
      // lineal con plasticidad
      material_list[i].E  = params[i*MAX_MAT_PARAM + 0];
      material_list[i].nu = params[i*MAX_MAT_PARAM + 1];
      material_list[i].plasticity = true;
      material_list[i].damage     = false;
    } else if (mat_types[i] == 2) {
      // lineal con daño
      material_list[i].E  = params[i*MAX_MAT_PARAM + 0];
      material_list[i].nu = params[i*MAX_MAT_PARAM + 1];
      material_list[i].plasticity = false;
      material_list[i].damage     = true;
    }
  }

  solver.max_its = 1000;
  solver.min_tol = 1.0e-8;
  this->dim = dim;
  if (dim == 2) {
    nx = size[0];
    ny = size[1];
    nvoi = 3;
    nn = nx * ny;
    npe = 4;
  } else if (dim == 3) {
    nx = size[0];
    ny = size[1];
    nz = size[2];
    nvoi = 6;
    nn = nx * ny * nz;
    npe = 8;
  }

  dx = lx/(nx-1);
  dy = ly/(ny-1);
  if (dim == 2) {
    nelem = (nx-1) * (ny-1);
  } else if (dim == 3) {
    dz = lz/(nz-1);
    nelem = (nx-1) * (ny-1) * (nz-1);
  }

  b  = (double*)malloc(nn*dim*sizeof(double));
  du = (double*)malloc(nn*dim*sizeof(double));
  u  = (double*)malloc(nn*dim*sizeof(double));
  stress = (double*)malloc(nelem*nvoi*sizeof(double));
  strain = (double*)malloc(nelem*nvoi*sizeof(double));
  elem_type = (int*)malloc(nelem*sizeof(int));


  for (int i=0; i<nn*dim; i++)
    u[i] = 0.0;

  if (dim == 2) {
    for (int e=0; e<nelem; e++)
      elem_type[e] = getElemType(e);
  } else {
    for (int ex=0; ex<nx-1; ex++) {
      for (int ey=0; ey<ny-1; ey++) {
	for (int ez=0; ez<nz-1; ez++) {
	  int e = glo_elem3D(ex,ey,ez);
	  elem_type[e] = getElemType(ex,ey,ez);
	}
      }
    }
  }

  NewRap_Its = 3;
  NewRap_Tol = 1.0e-5;

  if (dim == 2) {
    ell_init_2D (A, dim, nx, ny);
  } else if (dim == 3) {
    ell_init_3D (A, dim, nx, ny, nz);
  }

  return;
}

Problem::~Problem (void)
{
  ell_free (A);
  free(b);
  free(du);
  free(u);
  free(stress);
  free(strain);
  free(elem_type);
}

int Problem::getElemType (int e)
{
  if (distance(e) < 0.2) {
    return 1;
  }
  else {
    return 0;
  }
  return -1;

}

int Problem::getElemType (int ex, int ey, int ez)
{
  if (micro_type == 0) {
    // esfera en matriz
    double x1 = ex*dx + dx/2;
    double y1 = ey*dx + dy/2;
    double z1 = ez*dx + dz/2;
    double x2 = lx/2;
    double y2 = ly/2;
    double z2 = lz/2;
    double rad = micro_params[3];
    if ( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) < rad*rad )
      return 1;
    else 
      return 0;

  } else if (micro_type == 1) {

  }
}

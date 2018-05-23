#include <vector>
#include <iostream>
#include "micro.h"

Problem::Problem (int dim, int size[3], int micro_type, double *micro_params, int *mat_types, double *params)
{
  this->micro_type = micro_type;

  int nParams;
  if (micro_type == 0) {
    // mat 1 = matrix
    // mat 2 = sphere
    numMaterials = 2;
    nParams = 4;
  } else if (micro_type == 1) {
    // mat 1 = layer 1
    // mat 2 = layer 2
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
    material_list[i].E  = params[i*MAX_MAT_PARAM + 0];
    material_list[i].nu = params[i*MAX_MAT_PARAM + 1];
    if (mat_types[i] == 0) {
      // lineal
      material_list[i].plasticity = false;
      material_list[i].damage     = false;
    } else if (mat_types[i] == 1) {
      // con plasticidad
      material_list[i].plasticity = true;
      material_list[i].damage     = false;
    } else if (mat_types[i] == 2) {
      // con daño
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
    for (int ex=0; ex<nx-1; ex++) {
      for (int ey=0; ey<ny-1; ey++) {
	int e = glo_elem3D(ex,ey,0);
	elem_type[e] = getElemType(ex,ey);
      }
    }
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

int Problem::getElemType (int ex, int ey)
{
  if (micro_type == 0) {
    // esfera en matriz
    double x1 = ex*dx + dx/2;
    double y1 = ey*dy + dy/2;
    double x2 = lx/2;
    double y2 = ly/2;
    double rad = micro_params[3];
    if ( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) < rad*rad )
      return 1;
    else 
      return 0;

  } else if (micro_type == 1) {
    double y = ey*dy + dy/2;
    double espesor = micro_params[3];
    if (y < espesor)
      return 1;
    else 
      return 0;
  }
}

int Problem::getElemType (int ex, int ey, int ez)
{
  if (micro_type == 0) {
    // esfera en matriz
    double x1 = ex*dx + dx/2;
    double y1 = ey*dy + dy/2;
    double z1 = ez*dz + dz/2;
    double x2 = lx/2;
    double y2 = ly/2;
    double z2 = lz/2;
    double rad = micro_params[3];
    if ( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) < rad*rad )
      return 1;
    else 
      return 0;

  } else if (micro_type == 1) {
    double y = ey*dy + dy/2;
    double espesor = micro_params[3];
    if (y < espesor)
      return 1;
    else 
      return 0;
  }
}

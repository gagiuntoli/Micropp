#include <stdlib.h> 
#include <iostream>
#include "micro.h"

static Problem* micro = NULL;

extern "C"
{
  void
  micro_create_(int* dim)
  {

//    int dim = 0 ; 
    int nx = 0 ; 
    int ny = 0; 
    int nz = 0; 
    int cg_its = 0; 
    double cg_tol = 0.0; 

    int size[3];
    size[0] = 0.0;
    size[1] = 0.0;
    size[2] = 0.0;

std::cout<< dim[0] << "\n"; 

    micro = new Problem(dim[0], size, cg_its, cg_tol); 

  }


} 

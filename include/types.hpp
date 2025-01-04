/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Jimmy Aguilar Mena <kratsbinovish@gmail.com>
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

#pragma once

#include <iostream>
#include <map>

typedef struct {
  /* Results from newton-raphson loop */
  int its = 0;
  int solver_its = 0;
  bool converged = false;

  void print() {
    cout << "newton.its        : " << its << endl;
    cout << "newton.solver_its : " << solver_its << endl;
    cout << "newton.converged  : " << converged << endl;
  }

} newton_t;

typedef struct {
  int ngp = 1;
  int size[3];
  int type = 0;
  double geo_params[4] = {0.1, 0.1, 0.1, 0.1};
  struct material_base materials[4];
  int *coupling = nullptr;
  bool subiterations = false;
  int nsubiterations = 10;
  int mpi_rank = 0;
  int nr_max_its = NR_MAX_ITS;
  double nr_max_tol = NR_MAX_TOL;
  double nr_rel_tol = NR_REL_TOL;
  int cg_max_its = CG_MAX_ITS;
  double cg_abs_tol = CG_ABS_TOL;
  double cg_rel_tol = CG_REL_TOL;
  bool calc_ctan_lin = true;
  bool use_A0 = false;
  int its_with_A0 = 1;
  bool lin_stress = true;
  bool write_log = false;

  void print() {
    cout << "ngp  : " << ngp << endl;
    cout << "size : " << size[0] << endl;
    cout << "type  : " << type << endl;
    cout << "geo_params : " << geo_params[0] << endl;
    cout << "subiterations : " << subiterations << endl;
    cout << "nsubiterations : " << nsubiterations << endl;
    cout << "mpi_rank : " << mpi_rank << endl;
    cout << "nr_max_its : " << nr_max_its << endl;
    cout << "nr_max_tol : " << nr_max_tol << endl;
    cout << "nr_rel_tol : " << nr_rel_tol << endl;
    cout << "calc_ctan_lin : " << calc_ctan_lin << endl;
    cout << "use_A0 : " << use_A0 << endl;
    cout << "its_with_A0 : " << its_with_A0 << endl;
    cout << "lin_stress : " << lin_stress << endl;
    cout << "write_log : " << write_log << endl;
  }

} micropp_params_t;

enum {
  MIC_HOMOGENEOUS,
  MIC_SPHERE,
  MIC_LAYER_Y,
  MIC_CILI_FIB_X,
  MIC_CILI_FIB_Z,
  MIC_CILI_FIB_XZ,
  MIC_QUAD_FIB_XYZ,
  MIC_QUAD_FIB_XZ,
  MIC_QUAD_FIB_XZ_BROKEN_X,
  MIC3D_SPHERES,
  MIC3D_8,
  MIC3D_FIBS_20_ORDER,
  MIC3D_FIBS_20_DISORDER
};

static map<int, std::string> micro_names = {{MIC_HOMOGENEOUS, "MIC_HOMOGENEOUS"},
                                            {MIC_SPHERE, "MIC_SPHERE"},
                                            {MIC_LAYER_Y, "MIC_LAYER_Y"},
                                            {MIC_CILI_FIB_X, "MIC_CILI_FIB_X"},
                                            {MIC_CILI_FIB_Z, "MIC_CILI_FIB_Z"},
                                            {MIC_CILI_FIB_XZ, "MIC_CILI_FIB_XZ"},
                                            {MIC_QUAD_FIB_XYZ, "MIC_QUAD_FIB_XYZ"},
                                            {MIC_QUAD_FIB_XZ, "MIC_QUAD_FIB_XZ"},
                                            {MIC_QUAD_FIB_XZ_BROKEN_X, "MIC_QUAD_FIB_XZ_BROKEN_X"},
                                            {MIC3D_SPHERES, "MIC3D_SPHERES"},
                                            {MIC3D_8, "MIC3D_8"},
                                            {MIC3D_FIBS_20_ORDER, "MIC3D_FIBS_20_ORDER"},
                                            {MIC3D_FIBS_20_DISORDER, "MIC3D_FIBS_20_DISORDER"}};

enum { FE_LINEAR, FE_ONE_WAY, FE_FULL, MIX_RULE_CHAMIS };

static map<int, int> gp_counter = {{FE_LINEAR, 0}, {FE_ONE_WAY, 0}, {FE_FULL, 0}, {MIX_RULE_CHAMIS, 0}};

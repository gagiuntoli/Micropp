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

#define MAX_DIM 3
#define NUM_VAR_GP 7  // eps_p_1 (6) , alpha_1 (1)
#define MAX_MATERIALS 3

#define FILTER_REL_TOL 1.0e-5

#define D_EPS_CTAN_AVE 1.0e-8

#define CONSTXG 0.577350269189626

#define NR_MAX_TOL 1.0e-10
#define NR_MAX_ITS 4
#define NR_REL_TOL 1.0e-3  // factor against first residual

#define glo_elem(ex, ey, ez) ((ez) * (nx - 1) * (ny - 1) + (ey) * (nx - 1) + (ex))
#define intvar_ix(e, gp, var) ((e) * npe * NUM_VAR_GP + (gp) * NUM_VAR_GP + (var))

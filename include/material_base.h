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

#ifndef MATERIAL_BASE_H
#define MATERIAL_BASE_H

enum { MATERIAL_ELASTIC = 0, MATERIAL_PLASTIC, MATERIAL_DAMAGE };

#define D_EPS_CTAN 1.0e-8
#define SQRT_2DIV3 0.816496581

#ifdef __cplusplus
extern "C" {
#include <cstdio>
#else
#include <stdbool.h>
#include <stdio.h>
#endif

struct material_base {
  double E, nu, Ka, Sy;
  double k, mu, lambda;
  double Xt;
  int type;
};

void material_set(struct material_base *self, const int _type, const double E, const double nu, const double Ka,
                  const double Sy, const double Xt);

void material_print(const struct material_base *self);

#ifdef __cplusplus
}
#endif
#endif

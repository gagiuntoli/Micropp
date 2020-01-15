/*
 *  This source code is part of Micropp: a Finite Element library
 *  to solve composite materials micro-scale problems.
 *
 *  Copyright (C) - 2018
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
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

	void material_set(struct material_base *self, const int _type,
			  const double E, const double nu, const double Ka,
			  const double Sy, const double Xt);

	void material_print(const struct material_base *self);


#ifdef __cplusplus
}
#endif
#endif

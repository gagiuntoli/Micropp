/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                        Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include "material_base.h"

#include <cmath>
#include <cstdio>

struct material_t : public material_base {

	material_t() {
		E = NAN;
		nu = NAN;
		Ka = NAN;
		Sy = NAN;
		k = NAN;
		mu = NAN;
		lambda = NAN;
		type = -1;
		plasticity = false;
		damage = false;
	}

	void set(double _E, double _nu, double _Ka, double _Sy, int _type)
	{
		material_set(this, _E, _nu, _Ka, _Sy, _type);
	}

	void print() const
	{
		material_print(this);
	}
};

class material_damage : public material_t {

	public:
		material_damage(double _E, double _nu, double _Xt) {
			type = 2;
			E = _E;
			nu = _nu;
			Xt = _Xt;
		};

		void get_stress(const double eps[6], double stress[6], const double *params)
		{
			for (int i = 0; i < 3; ++i)
				stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
					    + 2 * mu * eps[i];

			for (int i = 3; i < 6; ++i)
				stress[i] = mu * eps[i];
		}
};


#endif

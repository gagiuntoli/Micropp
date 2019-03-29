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
#include <iostream>

using namespace std;

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

	material_t(const material_t * material) {
		E = material->E;
		nu = material->nu;
		Ka = material->Ka;
		Sy = material->Sy;
	}

	material_t * make_material(double * params, int type);
	static material_t * make_material(material_t *material);

	virtual void get_stress(const double eps[6], double stress[6], const double *history_params) {};

	virtual void print_n(void) {};

	void set(double _E, double _nu, double _Ka, double _Sy, int _type)
	{
		material_set(this, _type, _E, _nu, _Ka, _Sy, 0);
	}

	void print() const
	{
		material_print(this);
	}
};

class material_elastic : public material_t {

	public:
		material_elastic(double _E, double _nu) {
			E = _E;
			nu = _nu;
			Ka = 0;
			Sy = 0;
			Xt = 0;
		};

		void get_stress(const double eps[6], double stress[6], const double *history_params)
		{
			/* Elastic Material Law*/
			for (int i = 0; i < 3; ++i)
				stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
					    + 2 * mu * eps[i];

			for (int i = 3; i < 6; ++i)
				stress[i] = mu * eps[i];
		}

		void print_n(void) {
			cout << "Type : Elastic" << endl;
			cout << "E = " << E << " nu = " << nu << endl;
		}
};

class material_plastic : public material_t {

	public:
		material_plastic(double _E, double _nu, double _Ka, double _Sy) {
			E = _E;
			nu = _nu;
			Ka = _Ka;
			Sy = _Sy;
			Xt = 0;
		};

		void get_stress(const double eps[6], double stress[6], const double *history_params)
		{
			/* Elastic Material Law*/
			for (int i = 0; i < 3; ++i)
				stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
					    + 2 * mu * eps[i];

			for (int i = 3; i < 6; ++i)
				stress[i] = mu * eps[i];
		}

		void print_n(void) {
			cout << "Type : Plastic" << endl;
			cout << "E = " << E << " nu = " << nu << endl;
		}
};

class material_damage : public material_t {

	public:
		material_damage(double _E, double _nu, double _Xt) {
			E = _E;
			nu = _nu;
			Ka = 0;
			Sy = 0;
			Xt = _Xt;
		};

		void get_stress(const double eps[6], double stress[6], const double *history_params)
		{
			/* does not use history_params */
			for (int i = 0; i < 3; ++i)
				stress[i] = lambda * (eps[0] + eps[1] + eps[2]) \
					    + 2 * mu * eps[i];

			for (int i = 3; i < 6; ++i)
				stress[i] = mu * eps[i];
		}

		void print_n(void)
		{
			cout << "Type : Damage" << endl;
			cout << "E = " << E << " nu = " << nu
				<< " Xt = " << Xt << endl;
		}
};



#endif

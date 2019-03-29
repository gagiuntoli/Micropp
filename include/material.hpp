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
#include <cstring>
#include <iostream>

#define SQRT_2DIV3 0.816496581

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
	static material_t *make_material(material_t material);

	virtual void get_stress(const double *eps, double *stress, const double *history_params) const {};
	virtual void get_ctan(const double *eps, double *ctan, const double *history_params) const {};
	virtual bool evolute(const double *eps, const double *vars_old, double *vars_new) const {};

	virtual void print_n() const { cout << "I am base" << endl; };

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
			k = _E / (3. * (1. - 2. * _nu));
			mu = _E / (2. * (1. + _nu));
			lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));
			Ka = 0;
			Sy = 0;
			Xt = 0;
		};

		void get_stress(const double *eps, double *stress, const double *history_params) const;
		void get_ctan(const double *eps, double *ctan, const double *history_params) const;
		bool evolute(const double *eps, const double *vars_old, double *vars_new) const;
		void print_n() const;

};

class material_plastic : public material_t {

	public:
		material_plastic(double _E, double _nu, double _Ka, double _Sy) {
			E = _E;
			nu = _nu;
			k = _E / (3. * (1. - 2. * _nu));
			mu = _E / (2. * (1. + _nu));
			lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));
			Ka = _Ka;
			Sy = _Sy;
			Xt = 0;
		};

		void get_stress(const double *eps, double *stress, const double *history_params) const;
		void get_ctan(const double *eps, double *ctan, const double *history_params) const;
		bool evolute(const double *eps, const double *vars_old, double *vars_new) const;
		void print_n() const;

	private:
		bool plastic_law(const double eps[6],
				 const double *_eps_p_old,
				 const double *_alpha_old,
				 double *_dl,
				 double _normal[6],
				 double _s_trial[6]) const;
};

class material_damage : public material_t {

	public:
		material_damage(double _E, double _nu, double _Xt) {
			E = _E;
			nu = _nu;
			k = _E / (3. * (1. - 2. * _nu));
			mu = _E / (2. * (1. + _nu));
			lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));
			Ka = 0;
			Sy = 0;
			Xt = _Xt;
		};

		void get_stress(const double *eps, double *stress, const double *history_params) const;
		void get_ctan(const double *eps, double *ctan, const double *history_params) const;
		bool evolute(const double *eps, const double *vars_old, double *vars_new) const;
		void print_n() const;

};


#endif

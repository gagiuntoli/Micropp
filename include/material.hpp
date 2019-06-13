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


using namespace std;


struct material_t : public material_base {

	static material_t *make_material(const struct material_base material);

	virtual void init_vars(double *vars_old) const = 0;

	virtual void get_stress(const double *eps, double *stress,
				const double *history_params) const = 0;

	virtual void get_ctan(const double *eps, double *ctan,
			      const double *history_params) const = 0;

	/*
	 * <evolute> function updates <vars_new> using <eps> and
	 * <vars_old> returns <false> if the materials remains in
	 * linear zone and <true> if it has entered in non-linear zone.
	 */
	virtual bool evolute(const double *eps, const double *vars_old,
			     double *vars_new) const = 0;

	virtual void print() const = 0;

	protected:
	/*
	 * Apply_perturbation calculates <ctan> by applying the pertubation
	 * procedure.
	 */
	void apply_perturbation(const double *eps, double *ctan,
				const double *vars_old) const;


};


class material_elastic : public material_t {

	public:
		material_elastic(double _E, double _nu) {
			E = _E;
			nu = _nu;
			k = _E / (3. * (1. - 2. * _nu));
			mu = _E / (2. * (1. + _nu));
			lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));
			Ka = -1.0;
			Sy = -1.0;
			Xt = -1.0;
		};

		void init_vars(double *vars_old) const;

		void get_stress(const double *eps, double *stress,
				const double *history_params) const;

		void get_ctan(const double *eps, double *ctan,
			      const double *history_params) const;

		bool evolute(const double *eps, const double *vars_old,
			     double *vars_new) const;

		void print() const;

};


class material_plastic : public material_t {

	public:
		material_plastic(double _E, double _nu, double _Ka, double _Sy)
		{
			E = _E;
			nu = _nu;
			k = _E / (3. * (1. - 2. * _nu));
			mu = _E / (2. * (1. + _nu));
			lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));
			Ka = _Ka;
			Sy = _Sy;
			Xt = -1.0;
		};

		void init_vars(double *vars_old) const;

		void get_stress(const double *eps, double *stress,
				const double *history_params) const;

		void get_ctan(const double *eps, double *ctan,
			      const double *history_params) const;

		bool evolute(const double *eps, const double *vars_old,
			     double *vars_new) const;
		void print() const;

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
		material_damage(double _E, double _nu, double _Xt)
		{
			E = _E;
			nu = _nu;
			k = _E / (3. * (1. - 2. * _nu));
			mu = _E / (2. * (1. + _nu));
			lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));
			Ka = -1.0;
			Sy = -1.0;
			Xt = _Xt;
		};

		void init_vars(double *vars_old) const;

		void get_stress(const double *eps, double *stress,
				const double *history_params) const;

		void get_ctan(const double *eps, double *ctan,
			      const double *history_params) const;

		bool evolute(const double *eps, const double *vars_old,
			     double *vars_new) const;

		void print() const;

	private:
		double hardening_law(const double r) const;
		bool damage_law(const double *eps, const double e_old,
				const double D_old, double *_e, double *_D,
				double *stress_lin) const;

};


#endif

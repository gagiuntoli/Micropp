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
 *
 *************************************************************
 *
 * This code was though only for GPUs do to the current incompatibility
 * between OpenACC and the "virtual functions".
 *
 */

#ifndef MATERIAL_ACC_HPP
#define MATERIAL_ACC_HPP


#include "material_base.h"


#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>


using namespace std;


struct material_acc : public material_base {

	material_acc(const material_base material);

	void get_stress(const double *eps, double *stress,
			const double *history_params) const;

	void get_ctan(const double *eps, double *ctan,
		      const double *history_params) const;

	/*
	 * <evolute> function updates <vars_new> using <eps> and
	 * <vars_old> returns <false> if the materials remains in
	 * linear zone and <true> if it has entered in non-linear zone.
	 */
	bool evolute(const double *eps, const double *vars_old,
		     double *vars_new) const;

	void print();

};

void get_stress_elastic(const material_acc *material,
			const double *eps, double *stress,
			const double *history_params);

void get_ctan_elastic(const material_acc *material,
		      const double *eps, double *ctan,
		      const double *history_params);

bool evolute_elastic(const double *eps, const double *vars_old,
		     double *vars_new);

void get_stress_plastic(const material_acc *material,
			const double *eps, double *stress,
			const double *history_params);

void get_ctan_plastic(const material_acc *material,
		      const double *eps, double *ctan,
		      const double *history_params);

bool evolute_plastic(const material_acc *material, const double *eps, const double *vars_old,
		     double *vars_new);

void get_stress_damage(const material_acc *material,
		       const double *eps, double *stress,
		       const double *history_params);

void get_ctan_damage(const material_acc *material,
		     const double *eps, double *ctan,
		     const double *history_params);

bool evolute_damage(const material_acc *material, const double *eps, const double *vars_old,
		    double *vars_new);

#endif

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

#include <cmath>
#include <cstdio>

struct material_t {

	double E, nu, Ka, Sy;
	double k, mu, lambda;
	int type;
	bool plasticity, damage;

	material_t() :
		E(NAN), nu(NAN), Sy(NAN), Ka(NAN),
		k(NAN), mu(NAN), lambda(NAN),
		type(-1),
		plasticity(false), damage(false){}

	void set(double _E, double _nu, double _Ka, double _Sy, int _type)
	{
		E = _E;
		nu = _nu;
		Ka = _Ka;
		Sy = _Sy;
		type = _type;

		k = E / (3. * (1. - 2. * nu));
		mu = E / (2. * (1. + nu));
		lambda = nu * E / ((1. + nu) * (1. - 2. * nu));

		if (type == 0) {        // lineal
			plasticity = false;
			damage = false;
		} else if (type == 1) { // con plasticidad
			plasticity = true;
			damage = false;
		} else if (type == 2) {	// con daño
			plasticity = false;
			damage = true;
		}
	}

	void print()
	{
		printf("E = %lf, nu = %lf,	Ka = %lf, Sy = %lf, type = %d\n",
		       E, nu , Ka, Sy, type);
	}
};


#endif

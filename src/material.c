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


#include "material_base.h"


void material_set(struct material_base *self, const int _type,
		  const double _E, const double _nu, const double _Ka,
		  const double _Sy, const double _Xt)
{
	self->E = _E;
	self->nu = _nu;
	self->Ka = _Ka;
	self->Sy = _Sy;
	self->Xt = _Xt;
	self->type = _type;

	self->k = _E / (3. * (1. - 2. * _nu));
	self->mu = _E / (2. * (1. + _nu));
	self->lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));
}


void material_print(const struct material_base *self)
{
	printf("E = %e\nnu = %e\nKa = %e\nSy = %e\nXt = %e\nk = %e\nmu = %e\n"
	       "lambda = %e\ntype = %1d\n", self->E, self->nu , self->Ka,
	       self->Sy, self->Xt, self->k, self->mu, self->lambda, self->type);
}

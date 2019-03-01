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

#ifndef NEWTON_HPP
#define NEWTON_HPP


#define NR_MAX_TOL      1.0e-10
#define NR_MAX_ITS      4
#define NR_REL_TOL      1.0e-3 // factor against first residual


typedef struct {

	/* For getting performance results from newton-raphson loop */
	int its;
	double norms[NR_MAX_ITS] = { 0.0 };
	int solver_its[NR_MAX_ITS] = { 0 };
	double solver_norms[NR_MAX_ITS] = { 0.0 };

} newton_t;


#endif

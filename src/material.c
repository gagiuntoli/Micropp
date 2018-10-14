/*
 * This file is part of the Fortran_C++ program.
 * Copyright (c) 2018 Jimmy Aguilar Mena.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

typedef struct material_t {
    double E, nu, Sy, Ka;
    double k, mu, lambda;
    int type;
	bool plasticity, damage;
}

void set(material_t *material, double E, double nu, double Sy, double Ka, int type)
{
    material_set_(material, E, nu,  Sy, Ka, type);
}

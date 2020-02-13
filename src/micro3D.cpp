/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
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


#include "micropp.hpp"


template <>
void micropp<3>::set_displ_bc(const double eps[nvoi], double *u)
{
	const double eps_t[dim][dim] = {
		{       eps[0], 0.5 * eps[3], 0.5 * eps[4] },
		{ 0.5 * eps[3],       eps[1], 0.5 * eps[5] },
		{ 0.5 * eps[4], 0.5 * eps[5],       eps[2] } };


	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, 0); // z = 0
			const double coor[3] = { i * dx, j * dy, 0 };
			mvp(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, nz - 1); // z = lz
			const double coor[3] = { i * dx, j * dy, lz };
			mvp(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 0; k < nz; ++k) {
			const int n = nod_index3D(i, 0, k); // y = 0
			const double coor[3] = { i * dx, 0, k * dz };
			mvp(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 0; k < nz; ++k) {
			const int n = nod_index3D(i, ny - 1, k); // y = ly
			const double coor[3] = { i * dx, ly , k * dz };
			mvp(eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 0; j < ny; ++j) {
		for (int k = 0; k < nz; ++k) {
			const int n = nod_index3D(0, j, k); // x = 0
			const double coor[3] = { 0, j * dy , k * dz };
			mvp(eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 0; j < ny; ++j) {
		for (int k = 0; k < nz; ++k) {
			const int n = nod_index3D(nx - 1, j, k); // x = lx
			const double coor[3] = { lx, j * dy , k * dz };
			mvp(eps_t, coor, &u[n * dim]);
		}
	}
}


template <>
void micropp<3>::calc_bmat(int gp, double bmat[nvoi][npe * dim]) const
{
	const double dsh[8][3] = {
		{
			-(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			-(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			-(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			+(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			+(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			-(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			-(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			+(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			+(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz } };

	for (int i = 0; i < 8; ++i) {
		bmat[0][i * dim    ] = dsh[i][0];
		bmat[0][i * dim + 1] = 0;
		bmat[0][i * dim + 2] = 0;
		bmat[1][i * dim    ] = 0;
		bmat[1][i * dim + 1] = dsh[i][1];
		bmat[1][i * dim + 2] = 0;
		bmat[2][i * dim    ] = 0;
		bmat[2][i * dim + 1] = 0;
		bmat[2][i * dim + 2] = dsh[i][2];
		bmat[3][i * dim    ] = dsh[i][1];
		bmat[3][i * dim + 1] = dsh[i][0];
		bmat[3][i * dim + 2] = 0;
		bmat[4][i * dim    ] = dsh[i][2];
		bmat[4][i * dim + 1] = 0;
		bmat[4][i * dim + 2] = dsh[i][0];
		bmat[5][i * dim    ] = 0;
		bmat[5][i * dim + 1] = dsh[i][2];
		bmat[5][i * dim + 2] = dsh[i][1];
	}
}

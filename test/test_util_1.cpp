/*
 *  This is a test example for MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
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

#include <iostream>
#include <iomanip>

#include <ctime>
#include <cmath>
#include <cassert>

#include "util.hpp"

using namespace std;

int main (int argc, char *argv[])
{
	const double a1[2][2] = {
		{  10.66, -2.66 },
		{ -2.66,   10.66 } };
	const double v1[2] = { 47.1, -11.39 };
	const double s2_exact[2] = { 532.3834, -246.7034 };
	double s2[2];

	mvp_2(a1, v1, s2);
	for (int i = 0; i < 2; ++i)
		assert(fabs(s2[i] - s2_exact[i]) < 1.0e-10);

	const double a3[3][3] = {
		{  10.66, -2.66   , 8.3},
		{  10.66, -2.66   , 1.9},
		{ -2.66,  -10.66  , 7.2} };
	const double v3[3] = { -1.22, 7.1, -11.39 };
	const double s3_exact[3] = { -126.4282, -53.53220, -154.4488 };
	double s3[3];

	mvp_3(a3, v3, s3);
	for (int i = 0; i < 3; ++i)
		assert(fabs(s3[i] - s3_exact[i]) < 1.0e-10);

	return 0;
}

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


class gp_t {
	public:
		int nr_its[7];
		double *int_vars_n;
		double *u_n;
		double *int_vars_k;
		double *u_k;
		double macro_strain[6];
		double macro_stress[6];
		double macro_ctan[36];
		double nr_err[7];
		double inv_max;
		bool allocated;

		gp_t():
			int_vars_n(nullptr),
			u_n(nullptr),
			int_vars_k(nullptr),
			u_k(nullptr),
			inv_max(-1.0),
			allocated(false)
		{}

		~gp_t()
		{
			if (allocated) {
				free(int_vars_n);
				free(u_n);
				free(int_vars_k);
				free(u_k);
			}
		}

		void allocate(const int num_int_vars, const int nn, const int dim)
		{
			assert(!allocated);

			const int size = nn * dim;
			int_vars_n = (double *) calloc(num_int_vars, sizeof(double));
			int_vars_k = (double *) malloc(num_int_vars * sizeof(double));
			u_n = (double *) malloc(size * sizeof(double));
			u_k = (double *) malloc(size * sizeof(double));

			allocated = (int_vars_n && int_vars_k && u_n && u_k);
			assert(allocated);
		}


		void update_vars()
		{
			double *tmp = int_vars_n;
			int_vars_n = int_vars_k;
			int_vars_k = tmp;

			tmp = u_n;
			u_n = u_k;
			u_k = tmp;
		}
};

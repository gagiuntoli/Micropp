/*
 *  This source code is part of MicroPP: a finite element library
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
#include <string>
#include <sstream>
#include <cmath>
#include "micro.hpp"

using namespace std;

void micropp_t::output(int time_step, int gp_id)
{
	for (auto const &gp:gauss_list)
		if (gp.id == gp_id) {
			if (gp.int_vars_n != NULL)
				for (int i = 0; i < num_int_vars; ++i)
					vars_old[i] = gp.int_vars_n[i];
			else
				for (int i = 0; i < num_int_vars; ++i)
					vars_old[i] = 0.0;

			int nr_its;
			bool nl_flag;
			double nr_err;
			set_displ((double *)gp.macro_strain);
			newton_raphson(&nl_flag, &nr_its, &nr_err);

			calc_fields();
			write_vtu(time_step, gp_id);
			break;
		}
}

void micropp_t::write_vtu(int time_step, int gp_id)
{
	std::stringstream fname_vtu_s;
	fname_vtu_s << "micropp_" << gp_id << "_" << time_step << ".vtu";
	std::string fname_vtu = fname_vtu_s.str();

	ofstream file;
	file.open(fname_vtu);
	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"UnstructuredGrid\""
		 " version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	file << "<UnstructuredGrid>" << endl << "<Piece NumberOfPoints=\"" << nn;
	file << "\" NumberOfCells=\"" << nelem << "\">" << "\n<Points>" << endl;
	file << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\""
			 " format=\"ascii\">" << endl;

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				const double x = i * dx;
				const double y = j * dy;
				const double z = k * dz;
				file << x << " " << y << " " << z << endl;
			}
		}
	}
	file << "</DataArray>" << endl << "</Points>" << endl << "<Cells>" << endl;

	file << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int ez = 0; ez < nez; ez++) {
		for (int ey = 0; ey < ney; ey++) {
			for (int ex = 0; ex < nex; ex++) {
				const int n[] {
 					ez * (nx * ny) + ey * nx + ex,
					   ez * (nx * ny) + ey * nx + ex + 1,
					   ez * (nx * ny) + (ey + 1) * nx + ex + 1,
					   ez * (nx * ny) + (ey + 1) * nx + ex,
					   n[0] + (nx * ny),
					   n[1] + (nx * ny),
					   n[2] + (nx * ny),
					   n[3] + (nx * ny) } ;

				for (int i = 0; i < npe; ++i)
        			file << n[i] << ' ';
        		file << endl;
			}
		}
	}
	file << "</DataArray>" << endl;

	int ce = npe;
	file << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; e++) {
		file << ce << " ";
		ce += npe;
	}
	file << "\n</DataArray>" << endl;

	file << "<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	const int vtk_code = (dim == 2) ? 9 : 12;
	for (int e = 0; e < nelem; e++)
		file << vtk_code << " ";
	file << "\n</DataArray>" << endl;
	file << "</Cells>" << endl;

	file << "<PointData Vectors=\"displ,b\" >>" << endl;	// Vectors inside is a filter we should not use this here
	file << "<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
	for (int n = 0; n < nn; n++) {
		for (int d = 0; d < MAX_DIM; ++d)
			file << (dim == 2 && d == 2 ? 0.0 : u[n * dim + d]) << " ";
 		file << endl;
	}
	file << "</DataArray>" << endl;

	file << "<DataArray type=\"Float64\" Name=\"b\" NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
	for (int n = 0; n < nn; n++) {
		for (int d = 0; d < MAX_DIM; ++d)
			file << (dim == 2 && d == 2 ? 0.0 : b[n * dim + d]) << " ";
 		file << endl;
	}
	file << "</DataArray>" << endl;
	file << "</PointData>" << endl;

	file << "<CellData>" << endl;
	file << "<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; e++) {
		for (int v = 0; v < nvoi; v++)
			file << elem_strain[e * nvoi + v] << " ";
		file << endl;
	}
	file << "</DataArray>";

	file << "<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; e++) {
		for (int v = 0; v < nvoi; v++)
			file << elem_stress[e * nvoi + v] << " ";
		file << endl;
	}
	file << "</DataArray>";

	file << "<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; e++) {
		file << elem_type[e] << " ";
	}
	file << "\n</DataArray>" << endl;

	file << "<DataArray type=\"Float64\" Name=\"plasticity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; e++) {
		double plasticity = 0.0;
		if (vars_old != NULL)
			for (int gp = 0; gp < 8; gp++) {
				plasticity +=
					sqrt(vars_old[intvar_ix(e, gp, 0)] * vars_old[intvar_ix(e, gp, 0)] +
						    vars_old[intvar_ix(e, gp, 1)] * vars_old[intvar_ix(e, gp, 1)] +
						    vars_old[intvar_ix(e, gp, 2)] * vars_old[intvar_ix(e, gp, 2)] +
						    2 * vars_old[intvar_ix(e, gp, 3)] * vars_old[intvar_ix(e, gp, 3)] + 
							2 * vars_old[intvar_ix(e, gp, 4)] * vars_old[intvar_ix(e, gp, 4)] + 
							2 * vars_old[intvar_ix(e, gp, 5)] * vars_old[intvar_ix(e, gp, 5)]);
			}
		file << plasticity / 8 << " ";
	}
	file << "\n</DataArray>" << endl;

	file << "<DataArray type=\"Float64\" Name=\"hardening\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; e++) {
		double hardening = 0.0;
		if (vars_old != NULL)
			for (int gp = 0; gp < 8; gp++)
				hardening += vars_old[intvar_ix(e, gp, 6)];
		file << hardening / 8 << " ";
	}
	file << "\n</DataArray>" << endl;
	file << "</CellData>" << endl;
	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;

	file.close();
}

void micropp_t::write_info_files()
{
	ofstream file;
	if (output_files_header == false) {
		output_files_header = true;

		file.open("micropp_convergence.dat", std::ios_base::app);
		file << "# gp_id : ";
		for (auto const &gp:gauss_list)
			file << gp.id << " ";
		file << "\n# nl_flag [1] # inv_max [2] # inv_tol [3]" << endl << "# nr_its  [4] # nr_tol  [5]" << endl;
		file.close();

		file.open("micropp_eps_sig_ctan.dat", std::ios_base::app);
		file << "# gp_id : ";
		for (auto const &gp:gauss_list)
			file << gp.id << " ";
		file 
			<< "\n# epsxx [1] # epsyy [2] # epszz[3] # epsxy[4]  # epsxz[5]  # epsyz[6]"
			<< "\n# sigxx [7] # sigyy [8] # sigzz[9] # sigxy[10] # sigxz[11] # sigyz[12]" << endl;
		file.close();

		file.open("micropp_int_vars_n.dat", std::ios_base::app);
		file << "# gp_id : ";
		for (auto const &gp:gauss_list)
			file << gp.id << " ";
		file <<
 			"\n# eps_p_xx [1]  # eps_p_yy [2]  # eps_p_zz[3] "
 			"# eps_p_xy[4]  # eps_p_xz[5]  # eps_p_yz[6]  # alpha[7]"  << endl;
		file.close();
	}

	file.open("micropp_convergence.dat", std::ios_base::app);
	for (auto const &gp:gauss_list) {
		file << scientific;
		file << setw(3) << ((gp.int_vars_n == NULL) ? 0 : 1) << " ";
		file << setw(14) << gp.inv_max << " ";
		for (int i = 0; i < (1 + nvoi); ++i) {
			file << setw(14) << gp.nr_its[i] << " ";
			file << setw(14) << gp.nr_err[i] << " ";
		}
		file << " | ";
	}
	file << endl;
	file.close();

	file.open("micropp_eps_sig_ctan.dat", std::ios_base::app);
	for (auto const &gp:gauss_list) {
		for (int i = 0; i < 6; ++i)
			file << setw(14) << gp.macro_strain[i] << " ";
		for (int i = 0; i < 6; ++i)
			file << setw(14) << gp.macro_stress[i] << " ";
		for (int i = 0; i < 36; ++i)
			file << setw(14) << gp.macro_ctan[i] << " ";
		file << " | ";
	}
	file << endl;
	file.close();

	file.open("micropp_int_vars_n.dat", std::ios_base::app);
	for (auto const &gp:gauss_list) {
		for (int i = 0; i < num_int_vars; ++i)
			if (gp.int_vars_n != NULL)
				file << setw(14) << gp.int_vars_n[i] << " ";
			else
				file << setw(14) << 0.0 << " ";
		file << " | ";
	}
	file << endl;
	file.close();
}

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

template <int tdim>
void micropp<tdim>::output(int gp_id, const char *filename)
{
	INST_START;

	assert(gp_id < ngp);
	assert(gp_id >= 0);

	if (!gp_list[gp_id].allocated) {
		vars_old = vars_old_aux;
		vars_new = vars_new_aux;
		memset(vars_old, 0, num_int_vars * sizeof(double));
	} else {
		vars_old = gp_list[gp_id].int_vars_n;
		vars_new = gp_list[gp_id].int_vars_k;
	}

	double *u = gp_list[gp_id].u_k;

	calc_fields(u);
	write_vtu(u, filename);
}


template <int tdim>
void micropp<tdim>::write_vtu(const double *u, const char *filename)
{
	std::stringstream fname_vtu_s;
	fname_vtu_s << filename << ".vtu";
	std::string fname_vtu = fname_vtu_s.str();

	ofstream file;
	file.open(fname_vtu);

	file << "<?xml version=\"1.0\"?>\n"
		<<	"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
		<< "byte_order=\"LittleEndian\">\n"
		<< "<UnstructuredGrid>\n"
		<< "<Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\""
		<< nelem << "\">\n"
		<< "<Points>\n"
		<< "<DataArray type=\"Float64\" Name=\"Position\" "
		<< "NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				const double coor[3] = { i * dx, j * dy, k * dz };
				for (int d = 0; d < 3; d++)
					file << coor[d] << " ";
				file << endl;
			}
		}
	}
	file << "</DataArray>\n</Points>\n<Cells>\n" << endl;

	file
		<< "<DataArray type=\"Int32\" Name=\"connectivity\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {
				int n[8];
				get_elem_nodes(n, ex, ey, ez);
				for (int i = 0; i < npe; ++i)
					file << n[i] << ' ';
				file << endl;
			}
		}
	}
	file << "</DataArray>" << endl;

	int ce = npe;
	file
		<< "<DataArray type=\"Int32\" Name=\"offsets\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		file << ce << " ";
		ce += npe;
	}
	file << "\n</DataArray>" << endl;

	file
		<< "<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" "
		<< "format=\"ascii\">" << endl;
	const int vtk_code = (dim == 2) ? 9 : 12;
	for (int e = 0; e < nelem; ++e)
		file << vtk_code << " ";
	file << "\n</DataArray>" << endl;
	file << "</Cells>" << endl;

	file << "<PointData Vectors=\"displ\" >" << endl;
	file
		<< "<DataArray type=\"Float64\" Name=\"displ\" "
		<< "NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
	for (int n = 0; n < nn; ++n) {
		for (int d = 0; d < MAX_DIM; ++d)
			file << (dim == 2 && d == 2 ? 0.0 : \
				 ((isnan(u[n * dim + d]) ? -1 : u[n * dim + d]))) << " ";
		file << endl;
	}
	file << "</DataArray>" << endl;

	//  Deleted here
	file << "</PointData>" << endl;

	file << "<CellData>" << endl;
	file
		<< "<DataArray type=\"Float64\" Name=\"strain\" "
		<< "NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		for (int v = 0; v < nvoi; ++v)
			file << ((isnan(elem_strain[e * nvoi + v]) ? -1 : \
				  elem_strain[e * nvoi + v])) << " ";
		file << endl;
	}
	file << "</DataArray>\n";

	file
		<< "<DataArray type=\"Float64\" Name=\"stress\" "
		<< "NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		for (int v = 0; v < nvoi; ++v)
			file << ((isnan(elem_stress[e * nvoi + v]) ? -1 : \
				  elem_stress[e * nvoi + v])) << " ";
		file << endl;
	}
	file << "</DataArray>\n";

	file
		<< "<DataArray type=\"Int32\" Name=\"elem_type\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		file << elem_type[e] << " ";
	}
	file << "\n</DataArray>" << endl;

	file
		<< "<DataArray type=\"Float64\" Name=\"plasticity\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		double plasticity = 0.;
		for (int gp = 0; gp < npe; ++gp) {
			double tmp = 0.0;
			for (int v = 0; v < nvoi; ++v)
				tmp += vars_old[intvar_ix(e, gp, v)] *\
				       vars_old[intvar_ix(e, gp, v)];
			plasticity += sqrt(tmp);
		}
		file << plasticity / npe << " ";
	}
	file << "\n</DataArray>" << endl;

	file
		<< "<DataArray type=\"Float64\" Name=\"hardening\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		double hardening = 0.;
		for (int gp = 0; gp < npe; ++gp)
			hardening += vars_old[intvar_ix(e, gp, 6)];
		file << hardening / npe << " ";
	}
	file << "\n</DataArray>" << endl;
	file << "</CellData>" << endl;
	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;

	file.close();
}


template class micropp<2>;
template class micropp<3>;

/*
 *  This source code is part of MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
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

	calc_fields(gp_list[gp_id].u_k, gp_list[gp_id].vars_n);
	write_vtu(gp_list[gp_id].u_k, gp_list[gp_id].vars_n, filename);
}


template <int tdim>
void micropp<tdim>::output2(const int gp_id, const int elem_global,
			    const int time_step)
{
	/*
	 * This function writes the outfile
	 * "micropp-<elem_global>-<time_step>.vtu"
	 *
	 * gp_id       : is the local numeration of the Gauss point.
	 * elem_global : is the global numeration of th element in the
	 *               macro-scale.
	 * time_step   : is the time step in the macro-scale.
	 *
	 */
	INST_START;

	assert(gp_id < ngp);
	assert(gp_id >= 0);

	char filename[128];
	std::stringstream filename_stream;
	filename_stream << "micropp-" << elem_global << "-" << time_step;
	std::string file_name_string = filename_stream.str();
	strcpy(filename, file_name_string.c_str());

	calc_fields(gp_list[gp_id].u_k, gp_list[gp_id].vars_n);
	write_vtu(gp_list[gp_id].u_k, gp_list[gp_id].vars_n, filename);
}


template <int tdim>
void micropp<tdim>::write_vtu(double *u, double *vars_old, const char *filename)
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

	file << scientific;

	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				const double coor[3] = { i * dx, j * dy, k * dz };
				for (int d = 0; d < 3; d++)
					file << coor[d] << " ";
				file << "\n";
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
				file << "\n";
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
			file << ((dim == 2 && d == 2) ? 0.0 : u[n * dim + d]) << " ";
		file << "\n";
	}
	file << "</DataArray>" << endl;

	file << "</PointData>" << endl;

	file << "<CellData>" << endl;
	file
		<< "<DataArray type=\"Float64\" Name=\"strain\" "
		<< "NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		for (int v = 0; v < nvoi; ++v)
			file << elem_strain[e * nvoi + v] << " ";
		file << "\n";
	}
	file << "</DataArray>\n";

	file
		<< "<DataArray type=\"Float64\" Name=\"stress\" "
		<< "NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		for (int v = 0; v < nvoi; ++v)
			file << elem_stress[e * nvoi + v] << " ";
		file << "\n";
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
		double plasticity = 0.0;
		if (vars_old != nullptr) {
			for (int gp = 0; gp < npe; ++gp) {
				for (int v = 0; v < nvoi; ++v) {
					plasticity += vars_old[intvar_ix(e, gp, v)] *\
						      vars_old[intvar_ix(e, gp, v)];
				}
			}
		}
		plasticity += sqrt(plasticity);
		file << plasticity / npe << " ";
	}
	file << "\n</DataArray>" << endl;

	file
		<< "<DataArray type=\"Float64\" Name=\"damage_e\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		double damage = 0.0;
		if (vars_old != nullptr) {
			for (int gp = 0; gp < npe; ++gp) {
				damage += vars_old[intvar_ix(e, gp, 0)];
			}
		}
		file << damage / npe << " ";
	}
	file << "\n</DataArray>" << endl;

	file
		<< "<DataArray type=\"Float64\" Name=\"damage_D\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		double damage = 0.0;
		if (vars_old != nullptr) {
			for (int gp = 0; gp < npe; ++gp) {
				damage += vars_old[intvar_ix(e, gp, 1)];
			}
		}
		file << damage / npe << " ";
	}
	file << "\n</DataArray>" << endl;

	file
		<< "<DataArray type=\"Float64\" Name=\"hardening\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int e = 0; e < nelem; ++e) {
		double hardening = 0.;
		for (int gp = 0; gp < npe; ++gp) {
			if (vars_old != nullptr) {
				hardening += vars_old[intvar_ix(e, gp, 6)];
			} else {
				hardening += 0.0;
			}
		}
		file << hardening / npe << " ";
	}
	file << "\n</DataArray>" << endl;
	file << "</CellData>" << endl;
	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;

	file.close();
}


template <int tdim>
void micropp<tdim>::write_restart() const
{
	/*
	 *
	 * micropp-restart-<mpi_rank>-<#restart>.dat
	 */
	INST_START;

	int num_restart = 0;

	char filename[128];
	std::stringstream filename_stream;
	filename_stream << "micropp-restart-" << mpi_rank << "-" << num_restart
	         	<< ".bin";
	std::string file_name_string = filename_stream.str();
	strcpy(filename, file_name_string.c_str());

	ofstream file;
	//file.open (filename, ios::out | ios::binary);
	file.open (filename, ios::out);
	for (int igp = 0; igp < ngp; ++igp) {
		gp_list[igp].write_restart(file);
	}
	file.close();
}


template class micropp<3>;

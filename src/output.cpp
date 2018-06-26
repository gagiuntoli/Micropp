/*
 *  MicroPP : 
 *  Finite element library to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
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
#include "micro.h"

using namespace std;

void Problem::output (int time_step, int Gauss_ID, double *MacroStrain)
{
  list<GaussPoint_t>::iterator GaussPoint;
  for (GaussPoint=GaussPointList.begin(); GaussPoint !=  GaussPointList.end(); GaussPoint++) {
    if (GaussPoint->id == Gauss_ID) {
      if (GaussPoint->int_vars != NULL) {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = GaussPoint->int_vars[i];
      } else {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = 0.0;
      }
      break;
    }
  }

  double MacroStress[6];
  bool non_linear;
  setDisp(MacroStrain);
  newtonRaphson(&non_linear);

  calcDistributions();
  writeVtu(time_step, Gauss_ID);
}

void Problem::writeVtu (int time_step, int Gauss_ID)
{
  std::stringstream fname_vtu_s;
  fname_vtu_s  << "micropp_" << Gauss_ID << "_" << time_step << ".vtu";
  std::string fname_vtu  = fname_vtu_s.str();

  ofstream file;
  file.open (fname_vtu);
  file
      << "<?xml version=\"1.0\"?>" << endl
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl
      << "<UnstructuredGrid>"  << endl
      << "<Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << nelem << "\">" << endl
      << "<Points>" << endl
      << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

  double x, y, z;
  if (dim == 2) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
	x = i * dx;
	y = j * dy;
	file << x << " " << y << " "  << "0.0" << endl;
      }
    }
  }
  else if (dim == 3) {
    for (int k=0; k<nz; k++) {
      for (int j=0; j<ny; j++) {
	for (int i=0; i<nx; i++) {
	  x = i * dx;
	  y = j * dy;
	  z = k * dz;
	  file << x << " " << y << " "  << z << endl;
	}
      }
    }
  }
  file << "</DataArray>" << endl << "</Points>" << endl << "<Cells>" << endl;

  file << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  if (dim == 2) {
    for (int ey=0; ey<ny-1; ey++) {
      for (int ex=0; ex<nx-1; ex++) {
	int n0 = ey*nx     + ex;
	int n1 = ey*nx     + ex + 1;
	int n2 = (ey+1)*nx + ex + 1;
	int n3 = (ey+1)*nx + ex;
	file << n0 << " " << n1 << " " << n2 << " " << n3 << " " << endl;
      }
    }
  } else if (dim == 3) {
    for (int ez=0; ez<nz-1; ez++) {
      for (int ey=0; ey<ny-1; ey++) {
	for (int ex=0; ex<nx-1; ex++) {
	  int n0 = ez * (nx*ny) + ey*nx     + ex;
	  int n1 = ez * (nx*ny) + ey*nx     + ex + 1;
	  int n2 = ez * (nx*ny) + (ey+1)*nx + ex + 1;
	  int n3 = ez * (nx*ny) + (ey+1)*nx + ex;
	  int n4 = n0 + (nx*ny);
	  int n5 = n1 + (nx*ny);
	  int n6 = n2 + (nx*ny);
	  int n7 = n3 + (nx*ny);
	file
	  <<n0<<" " <<n1<< " " <<n2<< " " <<n3 << " "
	  <<n4<<" " <<n5<< " " <<n6<< " " <<n7 << " "<<endl;
	}
      }
    }
  }
  file << "</DataArray>" << endl;

  int ce = npe;
  file << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    file << ce << " ";
    ce += npe;
  }
  file << endl << "</DataArray>" << endl;

  file << "<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    if(dim == 2) {
      file << "9 ";
    } else if(dim == 3) {
      file << "12 ";
    }
  }
  file << endl;
  file << "</DataArray>" << endl << "</Cells>" << endl;

  file << "<PointData Vectors=\"displ,b\" >>" << endl; // Vectors inside is a filter we should not use this here
  file << "<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
  for (int n=0; n<nn; n++) {
    if (dim == 2) {
      file<<u[n*dim+0]<<" "<<u[n*dim + 1]<<" 0.0"<<endl;
    } else if (dim == 3) {
      file<<u[n*dim+0]<<" "<<u[n*dim + 1]<<" "<<u[n*dim + 2]<<endl;
    }
  }
  file << "</DataArray>" << endl;
  file << "<DataArray type=\"Float64\" Name=\"b\" NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
  for (int n=0; n<nn; n++) {
    if (dim == 2) {
      file<<b[n*dim+0]<<" "<<b[n*dim + 1]<<" 0.0"<<endl;
    } else if (dim == 3) {
      file<<b[n*dim+0]<<" "<<b[n*dim + 1]<<" "<<b[n*dim + 2]<<endl;
    }
  }
  file << "</DataArray>" << endl;
  file << "</PointData>" << endl;

  file << "<CellData>" << endl;

  file << "<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    for (int v=0; v<nvoi; v++)
      file << elem_strain[e*nvoi + v] << " ";
    file << endl;
  }
  file << "</DataArray>";

  file << "<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    for (int v=0; v<nvoi; v++)
      file << elem_stress[e*nvoi + v] << " ";
    file << endl;
  }
  file << "</DataArray>";

  file << "<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    file << elem_type[e] << " ";
  }
  file << endl << "</DataArray>" << endl;

  file << "<DataArray type=\"Float64\" Name=\"plasticity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    double plasticity = 0.0;
    if (vars_old != NULL) {
      if (dim == 3) {
	for (int gp=0; gp<8; gp++) {
	  plasticity += sqrt( \
	        vars_old[intvar_ix(e, gp, 0)] * vars_old[intvar_ix(e, gp, 0)]  + \
	        vars_old[intvar_ix(e, gp, 1)] * vars_old[intvar_ix(e, gp, 1)]  + \
	        vars_old[intvar_ix(e, gp, 2)] * vars_old[intvar_ix(e, gp, 2)]  + \
	      2*vars_old[intvar_ix(e, gp, 3)] * vars_old[intvar_ix(e, gp, 3)]  + \
	      2*vars_old[intvar_ix(e, gp, 4)] * vars_old[intvar_ix(e, gp, 4)]  + \
	      2*vars_old[intvar_ix(e, gp, 5)] * vars_old[intvar_ix(e, gp, 5)]);
	}
      }
    }
    file << plasticity/8 << " ";
  }
  file << endl << "</DataArray>" << endl;

  file << "<DataArray type=\"Float64\" Name=\"hardening\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    double hardening = 0.0;
    if (vars_old != NULL) {
      if (dim == 3) {
	for (int gp=0; gp<8; gp++) {
	  hardening += vars_old[intvar_ix(e, gp, 6)];
	}
      }
    }
    file << hardening/8 << " ";
  }
  file << endl << "</DataArray>" << endl;

  file << "</CellData>" << endl;
  file << "</Piece>" << endl << "</UnstructuredGrid>" << endl << "</VTKFile>" << endl;

  file.close();
}

void Problem::writeConvergenceFile (void)
{
  list <GaussPoint_t>::iterator GaussPoint;

  ofstream file;

  if (output_files_header == false) {
    output_files_header = true;

    file.open ("micropp_convergence.dat", std::ios_base::app);
    file << "GaussPointID : ";
    for (GaussPoint=GaussPointList.begin(); GaussPoint !=  GaussPointList.end(); GaussPoint++) {
      file << GaussPoint->id << " ";
    }
    file << endl;
    file.close();

    file.open ("micropp_eps_sig_ctan.dat", std::ios_base::app);
    file << "GaussPointID : ";
    for (GaussPoint=GaussPointList.begin(); GaussPoint !=  GaussPointList.end(); GaussPoint++) {
      file << GaussPoint->id << " ";
    }
    file << endl;
    file.close();
  }

  file.open ("micropp_convergence.dat", std::ios_base::app);
  for (GaussPoint=GaussPointList.begin(); GaussPoint !=  GaussPointList.end(); GaussPoint++) {
    file << scientific;
    file << setw(14) << GaussPoint->non_linear << " ";
    file << setw(14) << GaussPoint->convergence.I_reached << " ";
    file << setw(14) << GaussPoint->convergence.NR_Its_Stress << " ";
    file << setw(14) << GaussPoint->convergence.NR_Err_Stress << " ";
    for (int i=0; i<nvoi; i++) {
      file << setw(14) << GaussPoint->convergence.NR_Its_Ctan[i] << " ";
      file << setw(14) << GaussPoint->convergence.NR_Err_Ctan[i] << " ";
    }
    file << " | ";
  }
  file << endl;
  file.close();

  file.open ("micropp_eps_sig_ctan.dat", std::ios_base::app);
  for (GaussPoint=GaussPointList.begin(); GaussPoint !=  GaussPointList.end(); GaussPoint++) {
    for (int i=0; i<6; i++)
      file << setw(14) << GaussPoint->MacroStrain[i] << " ";
    for (int i=0; i<6; i++)
      file << setw(14) << GaussPoint->MacroStress[i] << " ";
    file << " | ";
  }
  file << endl;
  file.close();
}

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "micro.h"

using namespace std;

void Problem::output (int time_step, int elem, int macroGp_id, double *MacroStrain)
{
  double *int_vars;

  // search for the macro gauss point
  std::list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it !=  MacroGp_list.end(); it++) {
    if (it->id == macroGp_id) {
     int_vars = it->int_vars;
     break;
    }
  }
  if (it ==  MacroGp_list.end()) {
    cout << "output.cpp : Error the macro_id " << macroGp_id << " was not found in the list for plotting." << endl; 
  }

  double MacroStress[6];
  loc_hom_Stress(macroGp_id, MacroStrain, MacroStress);
  calcDistributions(int_vars);

  writeVtu(time_step, elem);
}

void Problem::writeVtu (int time_step, int elem)
{
  std::stringstream fname_pvtu_s, fname_vtu_s;;
  fname_pvtu_s << "micropp_" << elem << "_" << time_step << ".pvtu";
  fname_vtu_s  << "micropp_" << elem << "_" << time_step << ".vtu";
  std::string fname_pvtu = fname_pvtu_s.str();
  std::string fname_vtu  = fname_vtu_s.str();

  ofstream file;
  file.open (fname_pvtu);

  file 
    << "<?xml version=\"1.0\"?>" << endl
    << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl
    << "<PUnstructuredGrid GhostLevel=\"0\">" << endl
    << "<PPoints>"  << endl
    << "<PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\"/>"   << endl
    << "</PPoints>" << endl
    << "<PCells>"   << endl
    << "<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>" << endl
    << "<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>" << endl
    << "<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>" << endl
    << "</PCells>"  << endl
    << "<PPointData Vectors=\"displ\">" << endl
    << "<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />"  << endl
    << "<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />"  << endl
    << "</PPointData>" << endl
    << "<PCellData>"   << endl
    << "<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"" << 6 << "\"/>" << endl
    << "<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"" << 6 << "\"/>"    << endl
    << "<PDataArray type=\"Int32\"   Name=\"elem_type\" NumberOfComponents=\"1\"/>"  << endl
    << "</PCellData>"  << endl
    << "<Piece Source=\"" << fname_vtu << "\"/>" << endl
    << "</PUnstructuredGrid>" << endl
    << "</VTKFile>" << endl;
  file.close();

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
    } else if(dim == 3) {
      file<<u[n*dim+0]<<" "<<u[n*dim + 1]<<" "<<u[n*dim + 2]<<endl;
    }
  }
  file << "</DataArray>" << endl;
  file << "<DataArray type=\"Float64\" Name=\"b\" NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
  for (int n=0; n<nn; n++) {
    if (dim == 2) {
      file<<b[n*dim+0]<<" "<<b[n*dim + 1]<<" 0.0"<<endl;
    } else if(dim == 3) {
      file<<b[n*dim+0]<<" "<<b[n*dim + 1]<<" "<<b[n*dim + 2]<<endl;
    }
  }
  file << "</DataArray>" << endl;
  file << "</PointData>" << endl;

  file << "<CellData>" << endl;

  file << "<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    for (int v=0; v<nvoi; v++)
      file << strain[e*nvoi + v] << " ";
    file << endl;
  }
  file << "</DataArray>";

  file << "<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"" << nvoi << "\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    for (int v=0; v<nvoi; v++)
      file << stress[e*nvoi + v] << " ";
    file << endl;
  }
  file << "</DataArray>";

  file << "<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    file << elem_type[e] << " ";
  }
  file << endl << "</DataArray>" << endl;

  file << "</CellData>" << endl;
  file << "</Piece>" << endl << "</UnstructuredGrid>" << endl << "</VTKFile>" << endl;

  file.close();
}

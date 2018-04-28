#include <iostream>
#include <fstream>
#include <string>
#include "micro.h"

using namespace std;

void Problem::write_vtu (void)
{

  ofstream file;
  file.open ("solution.pvtu");

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
    << "<Piece Source=\"solution_0.vtu\"/>" << endl
    << "</PUnstructuredGrid>" << endl
    << "</VTKFile>" << endl;
  file.close();

  file.open ("solution_0.vtu");
  file
      << "<?xml version=\"1.0\"?>" << endl
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl
      << "<UnstructuredGrid>"  << endl
      << "<Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << nelem << "\">" << endl
      << "<Points>" << endl
      << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

  double x, y;
  for (int n=0; n<nn; n++) {
    if (dim == 2) {
      x = (n%nx) * dx;
      y = (n/nx) * dy;
      file << x << " " << y << " "  << "0.0" << endl;
    }
  }
  file << "</DataArray>" << endl << "</Points>" << endl << "<Cells>" << endl;

  file << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    if (dim == 2) {
      int xfactor = e%(nx-1);
      int yfactor = e/(ny-1);
      int n0 = yfactor     * nx + xfactor     ;
      int n1 = yfactor     * nx + xfactor + 1 ;
      int n2 = (yfactor+1) * nx + xfactor + 1 ;
      int n3 = (yfactor+1) * nx + xfactor     ;
      file << n0 << " " << n1 << " " << n2 << " " << n3 << " " << endl;
    }
  }
  file << "</DataArray>" << endl;

  int ce = npe;
  file << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e = 0 ; e < nelem ; e++) {
    file << ce << " ";
    ce += npe;
  }
  file << endl << "</DataArray>" << endl;

  file << "<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  for (int e=0; e<nelem; e++) {
    if(dim == 2) {
      file << "9 ";
    }
  }
  file << endl;
  file << "</DataArray>" << endl << "</Cells>" << endl;

  file << "<PointData Vectors=\"displ\">" << endl; // Vectors inside is a filter we should not use this here
  file << "<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >" << endl;
  for (int n=0; n<nn; n++) {
    if (dim == 2) {
      file << u[n*dim] << " " << u[n*dim + 1] << " 0.0" << endl;
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

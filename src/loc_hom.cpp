/*
 *  MicroPP : finite element library to solve microstructural problems for composite materials.
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

#include "micro.h"

void Problem::loc_hom_Ctan_Linear (double *Ctan)
{
  for (int i=0; i<nvoi; i++)
    for (int j=0; j<nvoi; j++)
      Ctan[i*nvoi + j] = CtanLinear[i][j];
}

void Problem::loc_hom_Stress_Linear (double *Strain, double *Stress)
{
  for (int i=0; i<nvoi; i++) {
    Stress[i] = 0.0;
    for (int j=0; j<nvoi; j++)
      Stress[i] += CtanLinear[i][j]*Strain[j];
  }
}

void Problem::loc_hom_Stress (int Gauss_ID, double *MacroStrain, double *MacroStress)
{
  bool not_allocated_yet = false;
  bool filter = true;
  double tol_filter = 1.0e-5;
  bool non_linear_history = false;

  // search for the macro gauss point, if not found just create and insert
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == Gauss_ID) {
      non_linear_history = it->non_linear;
      if (it->int_vars != NULL) {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = it->int_vars[i];
      } else {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = 0.0;
	not_allocated_yet = true;
      }
      for (int i=0; i<nvoi; i++)
	it->MacroStrain[i] = MacroStrain[i];
      break;
    }
  }

  if (it ==  MacroGp_list.end()) {
    MacroGp_t macroGp_new;
    macroGp_new.id = Gauss_ID;
    macroGp_new.non_linear = false;
    macroGp_new.non_linear_aux = false;
    macroGp_new.int_vars = NULL;
    macroGp_new.int_vars_aux = NULL;
    for (int i=0; i<nvoi; i++)
      for (int j=0; j<nvoi; j++)
	macroGp_new.CtanStatic[i*nvoi + j] = CtanLinear[i][j];
    MacroGp_list.push_back(macroGp_new);
    for (int i=0; i<num_int_vars; i++)
      vars_old[i] = 0.0;
    not_allocated_yet = true;
  }

  bool non_linear = false;

  if ((LinearCriteria(MacroStrain) == true) && (non_linear_history == false)) {

    for (int i=0; i<nvoi; i++) {
      MacroStress[i] = 0.0;
      for (int j=0; j<nvoi; j++)
	MacroStress[i] += CtanLinear[i][j] * MacroStrain[j];
      if (filter == true)
	if (fabs(MacroStress[i]) < tol_filter)
	  MacroStress[i] = 0.0;
    }
    for (int i=0; i<num_int_vars; i++)
      vars_new[i] = vars_old[i];

  } else {

    setDisp(MacroStrain);
    newtonRaphson(&non_linear);
    calcAverageStress(MacroStress);

    for (int v=0; v<nvoi; v++) {
      if (filter == true)
	if (fabs(MacroStress[v]) < tol_filter)
	  MacroStress[v] = 0.0;
    }
  }

  if (non_linear == true) {
    for (it=MacroGp_list.begin(); it !=  MacroGp_list.end(); it++) {
      if (it->id == Gauss_ID) {
	it->non_linear_aux = true;
	if (not_allocated_yet == true) {
	  it->int_vars = (double*)malloc(num_int_vars*sizeof(double));
	  it->int_vars_aux = (double*)malloc(num_int_vars*sizeof(double));
	}
	for (int i=0; i<num_int_vars; i++)
	  it->int_vars_aux[i] = vars_new[i];
	break;
      }
    }
  }
}

void Problem::loc_hom_Ctan (int Gauss_ID, double *MacroStrain, double *MacroCtan)
{
  bool filter = true;
  double tol_filter = 1.0e-2;
  bool non_linear_history = false;

  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == Gauss_ID) {
      non_linear_history = it->non_linear;
      if (it->int_vars != NULL) {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = it->int_vars[i];
      }
      break;
    }
  }
  if (it ==  MacroGp_list.end()) {
    for (int i=0; i<num_int_vars; i++)
      vars_old[i] = 0.0;
  }

  if (LinearCriteria(MacroStrain) == true) {

      for (int i=0; i<nvoi; i++)
	for (int j=0; j<nvoi; j++)
	  MacroCtan[i*nvoi + j] = CtanLinear[i][j];

  } else {

    double Strain_pert[6], Stress_0[6];
    double delta_Strain = 1.0e-10;
    bool non_linear;
    setDisp(MacroStrain);
    newtonRaphson(&non_linear);
    calcAverageStress(Stress_0);

    double stress_ave[6];
    for (int i=0; i<nvoi; i++) {
      for (int v=0; v<nvoi; v++)
	Strain_pert[v] = MacroStrain[v];
      Strain_pert[i] += delta_Strain;

      setDisp(Strain_pert);
      newtonRaphson(&non_linear);
      calcAverageStress(stress_ave);
      for (int v=0; v<nvoi; v++) {
	MacroCtan[v*nvoi + i] = (stress_ave[v] - Stress_0[v]) / delta_Strain;
	if (filter == true)
	  if (fabs(MacroCtan[v*nvoi + i]) < tol_filter)
	    MacroCtan[v*nvoi + i] = 0.0;
      }
    }
  }
}

void Problem::calcCtanLinear (void)
{
  double stress_ave[6], Strain_pert[6];
  double delta_Strain = 1.0e-8;
  bool non_linear;
  bool filter = true;
  double tol_filter = 1.0e-2;

  for (int i=0; i<nvoi; i++) {
    for (int v=0; v<nvoi; v++)
      Strain_pert[v] = 0.0;
    Strain_pert[i] += delta_Strain;

    setDisp(Strain_pert);
    newtonRaphson(&non_linear);
    calcAverageStress(stress_ave);

    for (int v=0; v<nvoi; v++) {
      CtanLinear[v][i] = stress_ave[v] / delta_Strain;
      if (filter == true)
	if (fabs(CtanLinear[v][i]) < tol_filter)
	  CtanLinear[v][i] = 0.0;
    }
  }
}

#define INV_MAX 10000

bool Problem::LinearCriteria (double *MacroStrain)
{
  double MacroStress[6];

  for (int i=0; i<nvoi; i++) {
    MacroStress[i] = 0.0;
    for (int j=0; j<nvoi; j++)
      MacroStress[i] += CtanLinear[i][j] * MacroStrain[j];
  }

  double Inv = Invariant_I1(MacroStress);
//  cout << " Inv = " << Inv << endl;
//  cout << " Stress = " 
//  << MacroStress[0] << " "
//  << MacroStress[1] << " " 
//  << MacroStress[2] << " "  
//  << MacroStress[3] << " " 
//  << MacroStress[4] << " " 
//  << MacroStress[5] << endl;

  if (fabs(Inv) > InvariantMax) InvariantMax = fabs(Inv);

  if (fabs(Inv) < INV_MAX) {
    LinCriteria = 1;
    return true;
  } else {
    LinCriteria = 0;
    return false;
  }
}

double Problem::Invariant_I1 (double *tensor)
{
  if (dim == 2)
    return tensor[0] + tensor[1];
  if (dim == 3)
    return tensor[0] + tensor[1] + tensor[2];
}

double Problem::Invariant_I2 (double *tensor)
{
  if (dim == 3)
    return tensor[0]*tensor[1] + tensor[0]*tensor[2] + tensor[1]*tensor[2] + tensor[3]*tensor[3] + tensor[4]*tensor[4] + tensor[5]*tensor[5];
}

void Problem::updateCtanStatic (void)
{
  bool filter = true;
  double tol_filter = 1.0e-2;

  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {

    if (LinearCriteria(it->MacroStrain) == true) {

      for (int i=0; i<nvoi; i++)
	for (int j=0; j<nvoi; j++)
	  it->CtanStatic[i*nvoi + j] = CtanLinear[i][j];

    } else {

      double Strain_pert[6], Stress_0[6];
      double delta_Strain = 1.0e-10;
      bool non_linear;
      setDisp(it->MacroStrain);
      newtonRaphson(&non_linear);
      calcAverageStress(Stress_0);

      double stress_ave[6];
      for (int i=0; i<nvoi; i++) {
	for (int v=0; v<nvoi; v++)
	  Strain_pert[v] = it->MacroStrain[v];
	Strain_pert[i] += delta_Strain;

	setDisp(Strain_pert);
	newtonRaphson(&non_linear);
	calcAverageStress(stress_ave);
	for (int v=0; v<nvoi; v++) {
	  it->CtanStatic[v*nvoi + i] = (stress_ave[v] - Stress_0[v]) / delta_Strain;
	  if (filter == true)
	    if (fabs(it->CtanStatic[v*nvoi + i]) < tol_filter)
	      it->CtanStatic[v*nvoi + i] = 0.0;
	}
      }
    }
    break;
  }
}

void Problem::getCtanStatic (int Gauss_ID, double *Ctan)
{
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == Gauss_ID) {
      for (int i=0; i<nvoi*nvoi; i++)
	Ctan[i] = it->CtanStatic[i];
      break;
    }
  }
}

void Problem::setMacroStrain(int Gauss_ID, double *MacroStrain)
{
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == Gauss_ID) {
      for (int i=0; i<nvoi; i++)
	it->MacroStrain[i] = MacroStrain[i];
      break;
    }
  }
  if (it ==  MacroGp_list.end()) {

    MacroGp_t macroGp_new;

    macroGp_new.id             = Gauss_ID;
    macroGp_new.non_linear     = false;
    macroGp_new.non_linear_aux = false;
    macroGp_new.int_vars       = NULL;
    macroGp_new.int_vars_aux   = NULL;

    MacroGp_list.push_back(macroGp_new);
  }
}

void Problem::getMacroStress(int Gauss_ID, double *MacroStress)
{
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == Gauss_ID) {
      for (int i=0; i<nvoi; i++)
	MacroStress[i] = it->MacroStress[i];
      break;
    }
  }
}

void Problem::getMacroCtan(int Gauss_ID, double *MacroCtan)
{
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == Gauss_ID) {
      for (int i=0; i<nvoi*nvoi; i++)
	MacroCtan[i] = it->MacroCtan[i];
      break;
    }
  }
}

void Problem::localizeHomogenize(void)
{
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->non_linear == false) {
      for (int i=0; i<num_int_vars; i++)
	vars_old[i] = 0.0;
    } else {
      for (int i=0; i<num_int_vars; i++)
	vars_old[i] = it->int_vars[i];
    }

    if ((LinearCriteria(it->MacroStrain) == true) && (it->non_linear == false)) {

      // S = CL : E
      for (int i=0; i<nvoi; i++) {
	it->MacroStress[i] = 0.0;
	for (int j=0; j<nvoi; j++)
	  it->MacroStress[i] += CtanLinear[i][j] * it->MacroStrain[j];
      }
      // C = CL
      for (int i=0; i<nvoi; i++) {
	for (int j=0; j<nvoi; j++)
	  it->MacroCtan[i*nvoi + j] = CtanLinear[i][j];
      }
      it->non_linear_aux = false;

    } else {

      // CALCULATE S (FEM)
      // HERE IS ONLY OPORTUNITY TO DETECT THE NON_LINEAR
      setDisp(it->MacroStrain);
      newtonRaphson(&it->non_linear_aux); 
      calcAverageStress(it->MacroStress);

      // CALCULATE C (FEM)
      // HERE WE DONT DETECT NON LINEARITY
      double Strain_1[6], Stress_0[6], Stress_1[6], dEps = 1.0e-10;
      bool non_linear_dummy;
      setDisp(it->MacroStrain);
      newtonRaphson(&non_linear_dummy);
      calcAverageStress(Stress_0);

      for (int i=0; i<nvoi; i++) {
	for (int v=0; v<nvoi; v++)
	  Strain_1[v] = it->MacroStrain[v];
	Strain_1[i] += dEps;

	setDisp(Strain_1);
	newtonRaphson(&non_linear_dummy);
	calcAverageStress(Stress_1);
	for (int v=0; v<nvoi; v++)
	  it->MacroCtan[v*nvoi + i] = (Stress_1[v] - Stress_0[v]) / dEps;
      }

    }
  }
}

void Problem::updateInternalVariables (void)
{
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {

    if (it->non_linear_aux == true) {

      if(it->non_linear == false) {
	it->non_linear = true;
	it->int_vars = (double *) malloc (num_int_vars * sizeof(double));
      }

      // WE ONLY NEED TO FILL <vars_new>
      bool non_linear_dummy;
      setDisp(it->MacroStrain);
      newtonRaphson(&non_linear_dummy);

      for (int i=0; i<num_int_vars; i++)
	it->int_vars[i] = vars_new[i];
    }
  }
}

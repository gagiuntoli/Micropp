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

#include "micro.h"

void Problem::calcCtanLinear()
{
  double Stress_1[6], Strain_1[6];
  double delta_Strain = 1.0e-8;
  bool non_linear;
  bool filter = true;
  double tol_filter = 1.0e-2;

  for (int i=0; i<nvoi; i++) {
    for (int v=0; v<nvoi; v++)
      Strain_1[v] = 0.0;
    Strain_1[i] += delta_Strain;

    setDisp(Strain_1);
    newtonRaphson(&non_linear);
    calcAverageStress(Stress_1);

    for (int v=0; v<nvoi; v++) {
      CtanLinear[v][i] = Stress_1[v] / delta_Strain;
      if ((filter == true) && (fabs(CtanLinear[v][i]) < tol_filter))
	CtanLinear[v][i] = 0.0;
    }
  }
}

bool Problem::LinearCriteria (double *MacroStrain)
{
  double MacroStress[6];

  for (int i=0; i<nvoi; i++) {
    MacroStress[i] = 0.0;
    for (int j=0; j<nvoi; j++) {
      MacroStress[i] += CtanLinear[i][j] * MacroStrain[j];
    }
  }
  double Inv = Invariant_I1(MacroStress);

  if (fabs(Inv) > I_reached) I_reached = fabs(Inv);

  return (fabs(Inv) < I_max) ? true : false;
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

void Problem::setMacroStrain(const int Gauss_ID, const double *MacroStrain)
{
  list<GaussPoint_t>::iterator GaussPoint;
  for (GaussPoint=GaussPointList.begin(); GaussPoint!=GaussPointList.end(); GaussPoint++) {
    if (GaussPoint->id == Gauss_ID) {
      for (int i=0; i<nvoi; i++)
	GaussPoint->MacroStrain[i] = MacroStrain[i];
      break;
    }
  }
  if (GaussPoint ==  GaussPointList.end()) {

    GaussPoint_t GaussPointNew;

    GaussPointNew.id             = Gauss_ID;
    GaussPointNew.non_linear     = false;
    GaussPointNew.non_linear_aux = false;
    GaussPointNew.int_vars       = NULL;
    GaussPointNew.int_vars_aux   = NULL;
    for (int i=0; i<nvoi; i++)
      GaussPointNew.MacroStrain[i] = MacroStrain[i];

    GaussPointNew.convergence.I_reached = -1.0e10;

    GaussPointList.push_back(GaussPointNew);
  }
}

void Problem::getMacroStress(const int Gauss_ID, double *MacroStress)
{
  list<GaussPoint_t>::iterator GaussPoint;
  for (GaussPoint=GaussPointList.begin(); GaussPoint!=GaussPointList.end(); GaussPoint++) {
    if (GaussPoint->id == Gauss_ID) {
      for (int i=0; i<nvoi; i++)
	MacroStress[i] = GaussPoint->MacroStress[i];
      break;
    }
  }
}

void Problem::getMacroCtan(const int Gauss_ID, double *MacroCtan)
{
  list<GaussPoint_t>::iterator GaussPoint;
  for (GaussPoint=GaussPointList.begin(); GaussPoint!=GaussPointList.end(); GaussPoint++) {
    if (GaussPoint->id == Gauss_ID) {
      for (int i=0; i<nvoi*nvoi; i++)
	MacroCtan[i] = GaussPoint->MacroCtan[i];
      break;
    }
  }
}

void Problem::localizeHomogenize()
{
  list<GaussPoint_t>::iterator GaussPoint;
  for (GaussPoint=GaussPointList.begin(); GaussPoint!=GaussPointList.end(); GaussPoint++) {
    if (GaussPoint->non_linear == false) {
      for (int i=0; i<num_int_vars; i++)
	vars_old[i] = 0.0;
    } else {
      for (int i=0; i<num_int_vars; i++)
	vars_old[i] = GaussPoint->int_vars[i];
    }

    I_reached = GaussPoint->convergence.I_reached;

    if ((LinearCriteria(GaussPoint->MacroStrain) == true) && (GaussPoint->non_linear == false)) {

      // S = CL : E
      for (int i=0; i<nvoi; i++) {
	GaussPoint->MacroStress[i] = 0.0;
	for (int j=0; j<nvoi; j++)
	  GaussPoint->MacroStress[i] += CtanLinear[i][j] * GaussPoint->MacroStrain[j];
      }
      // C = CL
      for (int i=0; i<nvoi; i++) {
	for (int j=0; j<nvoi; j++)
	  GaussPoint->MacroCtan[i*nvoi + j] = CtanLinear[i][j];
      }
      GaussPoint->non_linear_aux = false;
      GaussPoint->convergence.NR_Its_Stress = 0;
      GaussPoint->convergence.NR_Err_Stress = 0.0;
      for (int i=0; i<nvoi; i++) {
	GaussPoint->convergence.NR_Its_Ctan[i] = 0;
	GaussPoint->convergence.NR_Err_Ctan[i] = 0.0;
      }

    } else {

      // CALCULATE S (FEM)
      // HERE IS ONLY OPORTUNITY TO DETECT THE NON_LINEAR
      setDisp(GaussPoint->MacroStrain);
      newtonRaphson(&GaussPoint->non_linear_aux); 
      calcAverageStress(GaussPoint->MacroStress);

      GaussPoint->convergence.NR_Its_Stress = NR_its;
      GaussPoint->convergence.NR_Err_Stress = NR_norm;

      // CALCULATE C (FEM)
      // HERE WE DONT DETECT NON LINEARITY
      double Strain_1[6], Stress_0[6], Stress_1[6], dEps = 1.0e-10;
      bool non_linear_dummy;
      setDisp(GaussPoint->MacroStrain);
      newtonRaphson(&non_linear_dummy);
      calcAverageStress(Stress_0);

      for (int i=0; i<nvoi; i++) {
	for (int v=0; v<nvoi; v++)
	  Strain_1[v] = GaussPoint->MacroStrain[v];
	Strain_1[i] += dEps;

	setDisp(Strain_1);
	newtonRaphson(&non_linear_dummy);
	calcAverageStress(Stress_1);
	for (int v=0; v<nvoi; v++)
	  GaussPoint->MacroCtan[v*nvoi + i] = (Stress_1[v] - Stress_0[v]) / dEps;

	GaussPoint->convergence.NR_Its_Ctan[i] = NR_its;
	GaussPoint->convergence.NR_Err_Ctan[i] = NR_norm;
      }

    }
    
    GaussPoint->convergence.I_reached_aux = I_reached;
  }
}

void Problem::updateInternalVariables()
{
  list<GaussPoint_t>::iterator GaussPoint;
  for (GaussPoint=GaussPointList.begin(); GaussPoint!=GaussPointList.end(); GaussPoint++) {

    GaussPoint->convergence.I_reached = GaussPoint->convergence.I_reached_aux;

    if (GaussPoint->non_linear_aux == true) {

      if(GaussPoint->non_linear == false) {
	GaussPoint->non_linear = true;
	GaussPoint->int_vars = (double *) malloc (num_int_vars * sizeof(double));
      }

      // WE ONLY NEED TO FILL <vars_new>
      bool non_linear_dummy;
      setDisp(GaussPoint->MacroStrain);
      newtonRaphson(&non_linear_dummy);

      for (int i=0; i<num_int_vars; i++)
	GaussPoint->int_vars[i] = vars_new[i];
    }
  }
}

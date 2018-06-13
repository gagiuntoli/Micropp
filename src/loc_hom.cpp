#include "micro.h"

void Problem::loc_hom_Stress (int macroGp_id, double *MacroStrain, double *MacroStress)
{
  bool not_allocated_yet = false;
  bool filter = true;
  double tol_filter = 1.0e-5;
  bool non_linear_history = false;

  // search for the macro gauss point, if not found just create and insert
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == macroGp_id) {
      non_linear_history = it->non_linear;
      if (it->int_vars != NULL) {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = it->int_vars[i];
      } else {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = 0.0;
	not_allocated_yet = true;
      }
      break;
    }
  }

  if (it ==  MacroGp_list.end()) {
    MacroGp_t macroGp_new;
    macroGp_new.id = macroGp_id;
    macroGp_new.non_linear = false;
    macroGp_new.non_linear_aux = false;
    macroGp_new.int_vars = NULL;
    macroGp_new.int_vars_aux = NULL;
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
      if (it->id == macroGp_id) {
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

void Problem::updateIntVars (void)
{
  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    it->non_linear = it->non_linear_aux;
    if (it->int_vars != NULL) {
      for (int i=0; i<num_int_vars; i++)
	it->int_vars[i] = it->int_vars_aux[i];
    }
  }
}

void Problem::getNonLinearFlag (int macroGp_id, int *non_linear)
{
  *non_linear = 0;

  list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it!=MacroGp_list.end(); it++) {
    if (it->id == macroGp_id) {
      *non_linear = (it->non_linear == true) ? 1:0;
      break;
    }
  }
}

void Problem::loc_hom_Ctan (int macroGp_id, double *MacroStrain, double *MacroCtan)
{
  bool filter = true;
  double tol_filter = 1.0e-2;

  if (LinearCriteria(MacroStrain) == true) {

      for (int i=0; i<nvoi; i++)
	for (int j=0; j<nvoi; j++)
	  MacroCtan[i*nvoi + j] = CtanLinear[i][j];

  } else {

    list<MacroGp_t>::iterator it;
    for (it=MacroGp_list.begin(); it !=  MacroGp_list.end(); it++) {
      if (it->int_vars != NULL) {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = it->int_vars[i];
      } else {
	for (int i=0; i<num_int_vars; i++)
	  vars_old[i] = 0.0;
      }
    }

    if (it == MacroGp_list.end()) {
      for (int i=0; i<num_int_vars; i++)
	vars_old[i] = 0.0;
    }

    double Strain_pert[6], Stress_0[6];
    double delta_Strain = 0.00001;
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
  double delta_Strain = 0.00001;
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

#define I2_MAX 80000

bool Problem::LinearCriteria (double *MacroStrain)
{
  double MacroStress[6];

  for (int i=0; i<nvoi; i++) {
    MacroStress[i] = 0.0;
    for (int j=0; j<nvoi; j++)
      MacroStress[i] += CtanLinear[i][j] * MacroStrain[j];
  }

  double I2 = Invariant_I2(MacroStress);

  if (fabs(I2) < I2_MAX) {
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
  double J2;
  if (dim == 3)
    J2 = tensor[0]*tensor[1] + tensor[0]*tensor[2] + tensor[1]*tensor[2] + tensor[3]*tensor[3] + tensor[4]*tensor[4] + tensor[5]*tensor[5];
  return sqrt(J2);
}


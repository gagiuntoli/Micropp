#include "micro.h"

void Problem::loc_hom_Stress (int macroGp_id, double *MacroStrain, double *MacroStress)
{
  double *vars_old;
  bool allocate = false;

  // search for the macro gauss point, if not found just create and insert
  std::list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it !=  MacroGp_list.end(); it++) {
    if (it->id == macroGp_id) {
     cout << "Macro GP = "<< macroGp_id << " found. (loc_hom_Stress)" << endl; 
     if (it->int_vars != NULL) {
       vars_old = it->int_vars;
     }
     else {
       for (int i=0; i<(nelem*8*VARS_AT_GP); i++)
	 vars_dum_1[i] = 0.0;
       vars_old = vars_dum_1;
       allocate = true;
     }
     break;
    }
  }
  if (it ==  MacroGp_list.end()) {
    cout << "Macro GP = "<< macroGp_id << " NOT found, inserting GP... (loc_hom_Stress)" << endl; 
    MacroGp_t macroGp_new;
    macroGp_new.id = macroGp_id;
    macroGp_new.int_vars = NULL;
    MacroGp_list.push_back(macroGp_new);

    for (int i=0; i<(nelem*8*VARS_AT_GP); i++)
      vars_dum_1[i] = 0.0;

    vars_old = vars_dum_1;
    allocate = true;
  }

  double *vars_new = vars_dum_2;

  setDisp(MacroStrain);

  bool non_linear_flag = false;
  newtonRaphson(vars_old, vars_new, &non_linear_flag);

  if (non_linear_flag == true) {
    for (it=MacroGp_list.begin(); it !=  MacroGp_list.end(); it++) {
      if (it->id == macroGp_id) {

	if (allocate == true)
	  it->int_vars = (double*)malloc(nelem*8*VARS_AT_GP*sizeof(double));
	
	for (int i=0; i<(nelem*8*VARS_AT_GP); i++)
	  it->int_vars[i] = vars_new[i];

	break;
      }
    }
  }

  double stress_ave[6];
  calcAverageStress(vars_old, stress_ave);
  for (int v=0; v<nvoi; v++)
    MacroStress[v] = stress_ave[v];

}

void Problem::loc_hom_Ctan (int macroGp_id, double *MacroStrain, double *MacroCtan)
{
  double *vars_old;
  double *vars_new = vars_dum_2;

  // search for the macro gauss point, if not found just create and insert
  std::list<MacroGp_t>::iterator it;
  for (it=MacroGp_list.begin(); it !=  MacroGp_list.end(); it++) {
    if (it->id == macroGp_id) {
     vars_old = it->int_vars;
     break;
    }
  }
  if (it ==  MacroGp_list.end()) {

    for (int i=0; i<(nelem*8*VARS_AT_GP); i++)
      vars_dum_1[i] = 0.0;

    vars_old = vars_dum_1;
  }

  double Strain_pert[6], Stress_0[6];
  double delta_Strain = 0.00001;

  // calculate Stress_0 

  bool non_linear_flag;

  setDisp(MacroStrain);
  newtonRaphson(vars_old, vars_new, &non_linear_flag);

  double stress_ave[6];
  calcAverageStress(vars_old, stress_ave);
  for (int v=0; v<nvoi; v++)
    Stress_0[v] = stress_ave[v];

  // calculate Stress_1 .. 2 .. 3 , we calc directly MacroCtan

  for (int i=0; i<nvoi; i++) {

    for (int v=0; v<nvoi; v++)
      Strain_pert[v] = MacroStrain[v];

    Strain_pert[i] += delta_Strain; // we pertubate only one direction

    setDisp(Strain_pert);
    newtonRaphson(vars_old, vars_new, &non_linear_flag);

    calcAverageStress(vars_old, stress_ave);
    for (int v=0; v<nvoi; v++)
      MacroCtan[v*nvoi + i] = (stress_ave[v] - Stress_0[v]) / delta_Strain;

  }

}

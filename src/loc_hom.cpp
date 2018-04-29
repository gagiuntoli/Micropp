#include "micro.h"

void Problem::loc_hom_Stress (double *MacroStrain, double *MacroStress)
{
  setDisp(MacroStrain);
  newtonRaphson();
  calcAverageStress();
  for (int v=0; v<nvoi; v++)
    MacroStress[v] = stress_ave[v];
}

void Problem::loc_hom_Ctan (double *MacroStrain, double *MacroCtan)
{

  double Strain_pert[6], Stress_0[6];
  double delta_Strain = 0.00001;

  // calculate Stress_0 

  setDisp(MacroStrain);
  newtonRaphson();
  calcAverageStress();
  for (int v=0; v<nvoi; v++)
    Stress_0[v] = stress_ave[v];

  // calculate Stress_1 .. 2 .. 3 , we calc directly MacroCtan

  for (int i=0; i<nvoi; i++) {

    for (int v=0; v<nvoi; v++)
      Strain_pert[v] = MacroStrain[v];

    Strain_pert[i] += delta_Strain; // we pertubate only one direction

    setDisp(Strain_pert);
    newtonRaphson();
    calcAverageStress();
    for (int v=0; v<nvoi; v++)
      MacroCtan[v*nvoi + i] = (stress_ave[v] - Stress_0[v]) / delta_Strain;

  }

}

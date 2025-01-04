/*
 *  This source code is part of Micropp: a Finite Element library
 *  to solve composite materials micro-scale problems.
 *
 *  Copyright (C) - 2018
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

#include "material.hpp"

CUDA_HOSTDEV
material_t *material_t::make_material(const struct material_base material) {
  /*
   * This factory creates the corresponding subclass according to the
   * data members of <material>. This last is a "struct" than can be
   * manage also from C/Fortran easily and that is why we passed micropp
   * constructor the data in this form.
   */

  switch (material.type) {
    case 0:
      return new material_elastic(material.E, material.nu);
    case 1:
      return new material_plastic(material.E, material.nu, material.Ka, material.Sy);
    case 2:
      return new material_damage(material.E, material.nu, material.Xt);
    default:
      break;
  }
  return nullptr;
}

CUDA_HOSTDEV
void material_t::apply_perturbation(const double *eps, double *ctan, const double *vars_old) const {
  double stress_0[6];
  get_stress(eps, stress_0, vars_old);

  for (int i = 0; i < 6; ++i) {
    double eps_1[6];
    memcpy(eps_1, eps, 6 * sizeof(double));
    eps_1[i] += D_EPS_CTAN;

    double stress_1[6];
    get_stress(eps_1, stress_1, vars_old);

    for (int j = 0; j < 6; ++j) ctan[j * 6 + i] = (stress_1[j] - stress_0[j]) / D_EPS_CTAN;
  }
}

CUDA_HOSTDEV
void get_dev_tensor(const double tensor[6], double tensor_dev[6]) {
  memcpy(tensor_dev, tensor, 6 * sizeof(double));
  for (int i = 0; i < 3; i++) tensor_dev[i] -= (1 / 3.0) * (tensor[0] + tensor[1] + tensor[2]);
}

// ELASTIC MATERIAL

void material_elastic::init_vars(double *vars_old) const {}

CUDA_HOSTDEV
void material_elastic::get_stress(const double *eps, double *stress, const double *history_params) const {
  // stress[i][j] = lambda eps[k][k] * delta[i][j] + mu eps[i][j]
  for (int i = 0; i < 3; ++i) stress[i] = lambda * (eps[0] + eps[1] + eps[2]) + 2 * mu * eps[i];

  for (int i = 3; i < 6; ++i) stress[i] = mu * eps[i];
}

CUDA_HOSTDEV
void material_elastic::get_ctan(const double *eps, double *ctan, const double *history_params) const {
  // C = lambda * (1x1) + 2 mu I
  memset(ctan, 0, 6 * 6 * sizeof(double));

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) ctan[i * 6 + j] += lambda;

  for (int i = 0; i < 3; ++i) ctan[i * 6 + i] += 2 * mu;

  for (int i = 3; i < 6; ++i) ctan[i * 6 + i] = mu;
}

bool material_elastic::evolute(const double *eps, const double *vars_old, double *vars_new) const {
  // we don't have to evolute nothing is always linear
  return false;
}

void material_elastic::print() const {
  cout << "Type : Elastic" << endl;
  cout << scientific << "E = " << E << " nu = " << nu << endl;
}

// PLASTIC MATERIAL

void material_plastic::init_vars(double *vars_old) const {}

CUDA_HOSTDEV
bool material_plastic::plastic_law(const double eps[6], const double *_eps_p_old, const double *_alpha_old, double *_dl,
                                   double _normal[6], double _s_trial[6]) const {
  /*
   * Calculates _dl, _normal and _s_trial to used in the other plastic
   * material functions. Returns <true> if it enters in non-linear zone
   * <false> if not.
   */

  const double zeros[6] = {0.0};
  const double alpha_old = (_alpha_old) ? *_alpha_old : 0;
  const double *eps_p_old = (_eps_p_old) ? _eps_p_old : zeros;

  double eps_dev[6], eps_p_dev_1[6];

  get_dev_tensor(eps_p_old, eps_p_dev_1);
  get_dev_tensor(eps, eps_dev);

  for (int i = 0; i < 3; ++i) _s_trial[i] = 2 * mu * (eps_dev[i] - eps_p_dev_1[i]);

  for (int i = 3; i < 6; ++i) _s_trial[i] = mu * (eps_dev[i] - eps_p_dev_1[i]);

  double tmp = 0.0;
  for (int i = 0; i < 6; ++i) tmp += _s_trial[i] * _s_trial[i];
  double s_norm = sqrt(tmp);

  double f_trial = s_norm - SQRT_2DIV3 * (Sy + Ka * alpha_old);

  if (f_trial > 0) {
    for (int i = 0; i < 6; ++i) _normal[i] = _s_trial[i] / s_norm;
    *_dl = f_trial / (2. * mu * (1. + Ka / (3. * mu)));
    return true;

  } else {
    memset(_normal, 0, 6 * sizeof(double));
    *_dl = 0;
    return false;
  }
}

CUDA_HOSTDEV
void material_plastic::get_stress(const double *eps, double *stress, const double *history_params) const {
  double dl, normal[6], s_trial[6];
  const double *eps_p_old = (history_params) ? &(history_params[0]) : nullptr;
  const double *alpha_old = (history_params) ? &(history_params[6]) : nullptr;

  bool nl_flag = plastic_law(eps, eps_p_old, alpha_old, &dl, normal, s_trial);

  // sig_2 = s_trial + K * tr(eps) * 1 - 2 * mu * dl * normal;
  memcpy(stress, s_trial, 6 * sizeof(double));

  for (int i = 0; i < 3; ++i) stress[i] += k * (eps[0] + eps[1] + eps[2]);

  for (int i = 0; i < 6; ++i) stress[i] -= 2 * mu * dl * normal[i];
}

CUDA_HOSTDEV
void material_plastic::get_ctan(const double *eps, double *ctan, const double *vars_old) const {
  apply_perturbation(eps, ctan, vars_old);
}

bool material_plastic::evolute(const double *eps, const double *vars_old, double *vars_new) const {
  const double *eps_p_old = (vars_old) ? &(vars_old[0]) : nullptr;
  const double *alpha_old = (vars_old) ? &(vars_old[6]) : nullptr;
  double *eps_p_new = (vars_new) ? &(vars_new[0]) : nullptr;
  double *alpha_new = (vars_new) ? &(vars_new[6]) : nullptr;

  double dl, normal[6], s_trial[6];
  bool nl_flag = plastic_law(eps, eps_p_old, alpha_old, &dl, normal, s_trial);

  if (eps_p_old != nullptr && eps_p_new != nullptr)
    for (int i = 0; i < 6; ++i) eps_p_new[i] = eps_p_old[i] + dl * normal[i];

  if (alpha_old != nullptr && alpha_new != nullptr) *alpha_new = *alpha_old + SQRT_2DIV3 * dl + 0;

  return nl_flag;
}

void material_plastic::print() const {
  cout << "Type : Plastic" << endl;
  cout << "E = " << E << " nu = " << nu << " Ka = " << Ka << " Sy = " << Sy << endl;
}

/*
 * DAMAGE MATERIAL
 */

void material_damage::init_vars(double *vars_old) const {
  if (vars_old != nullptr) {
    const double Ey = 10.0e5;
    vars_old[0] = Ey / sqrt(E);
    vars_old[1] = 0.0;
  }
}

CUDA_HOSTDEV
double material_damage::hardening_law(const double r) const {
  const double Ey = 10.0e4;
  const double inf_Ey = 10. * Ey;

  const double H0 = 10.0;
  const double H1 = 5.0;

  const double r0 = Ey / sqrt(E);
  const double q0 = r0;                // strain_variable_init
  const double q1 = inf_Ey / sqrt(E);  // stress_variable_inf
  const double r1 = r0 + (q1 - q0) / H0;

  if (r < r0) {
    return 0.0;
  } else if (r >= r0 && r < r1) {
    return q0 + H0 * (r - r0);
  } else {
    return q1 + H1 * (r - r1);
  }
}

CUDA_HOSTDEV
bool material_damage::damage_law(const double *eps, const double r_old, const double D_old, double *_r, double *_D,
                                 double *_stress) const {
  /*
   * Calculates the linear stree <stress>, and <e> and <D> using the
   * strain <eps> and <e_old>
   *
   * e_old = vars_old[0]
   *
   * The function returns the real stress if stress != nullptr
   * if not returns D
   *
   */
  double stress_local[6] = {0.0};

  double *stress_ptr = (_stress != nullptr) ? _stress : stress_local;

  // First suppose we are in linear zone
  for (int i = 0; i < 3; ++i) stress_ptr[i] = lambda * (eps[0] + eps[1] + eps[2]) + 2 * mu * eps[i];

  for (int i = 3; i < 6; ++i) stress_ptr[i] = mu * eps[i];

  double product = 0.0;
  for (int i = 0; i < 6; ++i) {
    product += stress_ptr[i] * eps[i];
  }
  const double r = (product >= 0) ? sqrt(product) : 0;

  double r_old_a = (r_old < Xt / sqrt(E)) ? Xt / sqrt(E) : r_old;

  if (r <= r_old_a) {
    // Elastic
    *_r = r_old_a;
    *_D = D_old;
    return false;
  } else {
    // Inelastic (Damage)
    const double q = hardening_law(r);
    *_r = r;
    *_D = 1. - q / r;
    return true;
  }
}

CUDA_HOSTDEV
void material_damage::get_stress(const double *eps, double *stress, const double *vars_old) const {
  /*
   * Calculates the <stress> according to <eps> and <vars_old>.
   *
   * e_old = vars_old[0]
   *
   */

  const double r_old = (vars_old != nullptr) ? vars_old[0] : 0.0;
  const double D_old = (vars_old != nullptr) ? vars_old[1] : 0.0;
  double D, r;
  damage_law(eps, r_old, D_old, &r, &D, stress);

  for (int i = 0; i < 6; ++i) stress[i] *= (1 - D);
  // cout << D << endl;
}

CUDA_HOSTDEV
void material_damage::get_ctan(const double *eps, double *ctan, const double *vars_old) const {
  apply_perturbation(eps, ctan, vars_old);
}

bool material_damage::evolute(const double *eps, const double *vars_old, double *vars_new) const {
  /* Assign new values to <vars_new> according to <eps> and <vars_old>.
   * returns <true> if the material has entered in non-linear range,
   * <false> if not.
   */
  const double r_old = (vars_old) ? vars_old[0] : 0;
  const double D_old = (vars_old) ? vars_old[1] : 0;
  double *r_new = (vars_new) ? &(vars_new[0]) : nullptr;
  double *D_new = (vars_new) ? &(vars_new[1]) : nullptr;

  bool non_linear = damage_law(eps, r_old, D_old, r_new, D_new, nullptr);

  return non_linear;
}

void material_damage::print() const {
  cout << "Type : Damage" << endl;
  cout << "E = " << E << " nu = " << nu << " Xt = " << Xt << endl;
}

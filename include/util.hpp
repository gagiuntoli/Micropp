/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                        Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

inline uint64_t devest(const vector<uint64_t> &in, const uint64_t mean) {
  uint64_t out = 0;
  for (const auto &x : in) {
    const uint64_t tmp = (x - mean);
    out += tmp * tmp;
  }

  return sqrt(out / in.size());
}

constexpr int mypow(int v, int e) { return (e == 0) ? 1 : v * mypow(v, e - 1); }

inline void print_vec(const double *vec, int n, const char file_name[]) {
  FILE *file = fopen(file_name, "w");
  for (int i = 0; i < n; ++i) {
    fprintf(file, "[%lf]\n", vec[i]);
  }
  fclose(file);
}

template <typename T, int n>
inline void mvp(const T m[n][n], const T x[n], T *y) {
  for (int i = 0; i < n; ++i) {
    T tmp = 0.0;
    for (int j = 0; j < n; ++j) {
      tmp += m[i][j] * x[j];
    }
    y[i] = tmp;
  }
}

template <typename T, int n>
inline double norm(const T v1[n]) {
  /*
   * Returns sqrt (v1[0] * v1[0] + ... + v1[n-1] * v1[n-1])
   */

  T tmp = 0;
  for (int i = 0; i < n; ++i) {
    tmp += v1[i] * v1[i];
  }
  return sqrt((double)tmp);
}

template <typename T, int n>
inline T dot_prod(const T v1[n], const T v2[n]) {
  /*
   * Returns v1[0] * v2[0] + ... + v1[n-1] * v2[n-1]
   */

  T tmp = 0;
  for (int i = 0; i < n; ++i) {
    tmp += v1[i] * v2[i];
  }
  return tmp;
}

inline bool point_inside_sphere(const double center[3], const double radius, const double point[3]) {
  /*
   * Returns <true> if <point> is inside the sphere with <center> and
   * <radius>. Returns <false> otherwise.
   */

  const double v[3] = {point[0] - center[0], point[1] - center[1], point[2] - center[2]};

  return (norm<double, 3>(v) < radius) ? true : false;
}

inline bool point_inside_cilinder_inf(const double dir[3], const double center[3], const double radius,
                                      const double point[3]) {
  /*
   * Returns <true> if <point> is inside the infinite cilinder with
   * <direction>, <center> and <radius>. Returns <false> otherwise.
   */

  const double v[3] = {point[0] - center[0], point[1] - center[1], point[2] - center[2]};

  const double dir_dot_v = dot_prod<double, 3>(dir, v);
  const double norm_dir = norm<double, 3>(dir);
  const double norm_v = norm<double, 3>(v);
  const double cos_tetha = dir_dot_v / (norm_dir * norm_v);
  const double sin_tetha = sqrt(1 - cos_tetha * cos_tetha);
  const double d_1 = norm_v * sin_tetha;

  return (d_1 <= radius) ? true : false;
}

inline double invert_3x3(const double mat[3][3], double mat_inv[3][3]) {
  const double det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) -
                     mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]) +
                     mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);

  mat_inv[0][0] = +(mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) / det;
  mat_inv[0][1] = -(mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]) / det;
  mat_inv[0][2] = +(mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]) / det;

  mat_inv[1][0] = -(mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2]) / det;
  mat_inv[1][1] = +(mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2]) / det;
  mat_inv[1][2] = -(mat[0][0] * mat[2][1] - mat[2][0] * mat[0][1]) / det;

  mat_inv[2][0] = +(mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]) / det;
  mat_inv[2][1] = -(mat[0][0] * mat[1][2] - mat[1][0] * mat[0][2]) / det;
  mat_inv[2][2] = +(mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1]) / det;

  return det;
}

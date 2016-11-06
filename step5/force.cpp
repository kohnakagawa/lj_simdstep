#include <stdio.h>
#include <immintrin.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/stat.h>
//----------------------------------------------------------------------
const double density = 1.0;
//const double density = 0.5;
const int N = 400000;
const int MAX_PAIRS = 30 * N;
double L = 50.0;
const double dt = 0.001;
const int D = 4;
enum {X, Y, Z};
double q[N][D];
double p[N][D];

int particle_number = 0;
int number_of_pairs = 0;
int number_of_partners[N];
int i_particles[MAX_PAIRS];
int j_particles[MAX_PAIRS];
int pointer[N], pointer2[N];
int sorted_list[MAX_PAIRS];

const double CUTOFF_LENGTH = 3.0;
const double SEARCH_LENGTH = 3.3;
const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
//----------------------------------------------------------------------
typedef double v4df __attribute__((vector_size(32)));
//----------------------------------------------------------------------
void
print256(v4df r) {
  double *a = (double*)(&r);
  printf("%.10f %.10f %.10f %.10f\n", a[0], a[1], a[2], a[3]);
}
//----------------------------------------------------------------------
void
add_particle(double x, double y, double z) {
  static std::mt19937 mt(2);
  std::uniform_real_distribution<double> ud(0.0, 0.1);
  q[particle_number][X] = x + ud(mt);
  q[particle_number][Y] = y + ud(mt);
  q[particle_number][Z] = z + ud(mt);
  particle_number++;
}
//----------------------------------------------------------------------
double
myclock(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}
//----------------------------------------------------------------------
void
register_pair(int index1, int index2) {
  int i, j;
  if (index1 < index2) {
    i = index1;
    j = index2;
  } else {
    i = index2;
    j = index1;
  }
  i_particles[number_of_pairs] = i;
  j_particles[number_of_pairs] = j;
  number_of_partners[i]++;
  number_of_pairs++;
}
//----------------------------------------------------------------------
void
sortpair(void) {
  const int pn = particle_number;
  int pos = 0;
  pointer[0] = 0;
  for (int i = 0; i < pn - 1; i++) {
    pos += number_of_partners[i];
    pointer[i + 1] = pos;
  }
  for (int i = 0; i < pn; i++) {
    pointer2[i] = 0;
  }
  const int s = number_of_pairs;
  for (int k = 0; k < s; k++) {
    int i = i_particles[k];
    int j = j_particles[k];
    int index = pointer[i] + pointer2[i];
    sorted_list[index] = j;
    pointer2[i] ++;
  }
}
//----------------------------------------------------------------------
void
makepair(void) {
  const double SL2 = SEARCH_LENGTH * SEARCH_LENGTH;
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    number_of_partners[i] = 0;
  }
  for (int i = 0; i < particle_number - 1; i++) {
    for (int j = i + 1; j < particle_number; j++) {
      const double dx = q[i][X] - q[j][X];
      const double dy = q[i][Y] - q[j][Y];
      const double dz = q[i][Z] - q[j][Z];
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < SL2) {
        register_pair(i, j);
      }
    }
  }
}
//----------------------------------------------------------------------
void
init(void) {
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L / s);
  int sy = static_cast<int>(L / s);
  int sz = static_cast<int>(L / s);
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        double x = ix * s;
        double y = iy * s;
        double z = iz * s;
        add_particle(x     , y   , z);
        add_particle(x     , y + hs, z + hs);
        add_particle(x + hs  , y   , z + hs);
        add_particle(x + hs  , y + hs, z);
      }
    }
  }
  for (int i = 0; i < particle_number; i++) {
    p[i][X] = 0.0;
    p[i][Y] = 0.0;
    p[i][Z] = 0.0;
  }
}
//----------------------------------------------------------------------
void
force_pair(void) {
  for (int k = 0; k < number_of_pairs; k++) {
    const int i = i_particles[k];
    const int j = j_particles[k];
    double dx = q[j][X] - q[i][X];
    double dy = q[j][Y] - q[i][Y];
    double dz = q[j][Z] - q[i][Z];
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    p[i][X] += df * dx;
    p[i][Y] += df * dy;
    p[i][Z] += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
force_pair_swp(void) {
  int k = 0;
  int i_a = i_particles[k];
  int j_a = j_particles[k];
  double dx_b = q[j_a][X] - q[i_a][X];
  double dy_b = q[j_a][Y] - q[i_a][Y];
  double dz_b = q[j_a][Z] - q[i_a][Z];
  double dx_a, dy_a, dz_a;
  int i_b, j_b;
  double df;
  for (k = 1; k < number_of_pairs; k++) {
    dx_a = dx_b;
    dy_a = dy_b;
    dz_a = dz_b;
    i_b = i_particles[k];
    j_b = j_particles[k];
    dx_b = q[j_b][X] - q[i_b][X];
    dy_b = q[j_b][Y] - q[i_b][Y];
    dz_b = q[j_b][Z] - q[i_b][Z];
    const double r2 = (dx_a * dx_a + dy_a * dy_a + dz_a * dz_a);
    const double r6 = r2 * r2 * r2;
    df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    if (r2 > CL2) df = 0.0;
    p[i_a][X] += df * dx_a;
    p[i_a][Y] += df * dy_a;
    p[i_a][Z] += df * dz_a;
    p[j_a][X] -= df * dx_a;
    p[j_a][Y] -= df * dy_a;
    p[j_a][Z] -= df * dz_a;
    i_a = i_b;
    j_a = j_b;
  }
  dx_a = dx_b;
  dy_a = dy_b;
  dz_a = dz_b;
  const double r2 = (dx_a * dx_a + dy_a * dy_a + dz_a * dz_a);
  const double r6 = r2 * r2 * r2;
  df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
  if (r2 > CL2) df = 0.0;
  p[i_a][X] += df * dx_a;
  p[i_a][Y] += df * dy_a;
  p[i_a][Z] += df * dz_a;
  p[j_a][X] -= df * dx_a;
  p[j_a][Y] -= df * dy_a;
  p[j_a][Z] -= df * dz_a;
}
//----------------------------------------------------------------------
#define p4(x) printf("%.10f %.10f %.10f %.10f\n",dx_##x,dy_##x,dz_##x,0.0);
//----------------------------------------------------------------------
void
force_pair_intrin(void) {
  const v4df vzero = _mm256_set_pd(0, 0, 0, 0);
  const v4df vcl2 = _mm256_set_pd(CL2, CL2, CL2, CL2);
  const v4df vc24 = _mm256_set_pd(24 * dt, 24 * dt, 24 * dt, 24 * dt);
  const v4df vc48 = _mm256_set_pd(48 * dt, 48 * dt, 48 * dt, 48 * dt);
  int k = 0;
  int i_a1 = i_particles[k];
  int j_a1 = j_particles[k];
  int i_a2 = i_particles[k + 1];
  int j_a2 = j_particles[k + 1];
  int i_a3 = i_particles[k + 2];
  int j_a3 = j_particles[k + 2];
  int i_a4 = i_particles[k + 3];
  int j_a4 = j_particles[k + 3];
  v4df vqi_a1 = _mm256_load_pd((double*)(q + i_a1));
  v4df vqj_a1 = _mm256_load_pd((double*)(q + j_a1));
  v4df vdq_b1 = vqj_a1 - vqi_a1;

  v4df vqi_a2 = _mm256_load_pd((double*)(q + i_a2));
  v4df vqj_a2 = _mm256_load_pd((double*)(q + j_a2));
  v4df vdq_b2 = vqj_a2 - vqi_a2;

  v4df vqi_a3 = _mm256_load_pd((double*)(q + i_a3));
  v4df vqj_a3 = _mm256_load_pd((double*)(q + j_a3));
  v4df vdq_b3 = vqj_a3 - vqi_a3;

  v4df vqi_a4 = _mm256_load_pd((double*)(q + i_a4));
  v4df vqj_a4 = _mm256_load_pd((double*)(q + j_a4));
  v4df vdq_b4 = vqj_a4 - vqi_a4;

  v4df vdq_a1;
  v4df vdq_a2;
  v4df vdq_a3;
  v4df vdq_a4;

  int i_b1, j_b1;
  int i_b2, j_b2;
  int i_b3, j_b3;
  int i_b4, j_b4;
  v4df vdf;

  for (k = 4; k < (number_of_pairs) / 4 * 4; k += 4) {
    vdq_a1 = vdq_b1;
    vdq_a2 = vdq_b2;
    vdq_a3 = vdq_b3;
    vdq_a4 = vdq_b4;

    i_b1 = i_particles[k];
    j_b1 = j_particles[k];
    i_b2 = i_particles[k + 1];
    j_b2 = j_particles[k + 1];
    i_b3 = i_particles[k + 2];
    j_b3 = j_particles[k + 2];
    i_b4 = i_particles[k + 3];
    j_b4 = j_particles[k + 3];

    v4df vqi_b1 = _mm256_load_pd((double*)(q + i_b1));
    v4df vqj_b1 = _mm256_load_pd((double*)(q + j_b1));
    vdq_b1 = vqj_b1 - vqi_b1;

    v4df vqi_b2 = _mm256_load_pd((double*)(q + i_b2));
    v4df vqj_b2 = _mm256_load_pd((double*)(q + j_b2));
    vdq_b2 = vqj_b2 - vqi_b2;

    v4df vqi_b3 = _mm256_load_pd((double*)(q + i_b3));
    v4df vqj_b3 = _mm256_load_pd((double*)(q + j_b3));
    vdq_b3 = vqj_b3 - vqi_b3;

    v4df vqi_b4 = _mm256_load_pd((double*)(q + i_b4));
    v4df vqj_b4 = _mm256_load_pd((double*)(q + j_b4));
    vdq_b4 = vqj_b4 - vqi_b4;

    v4df vr2_1x = vdq_a1 * vdq_a1;
    v4df vr2_1y = _mm256_permute4x64_pd(vr2_1x, 201);
    v4df vr2_1z = _mm256_permute4x64_pd(vr2_1x, 210);
    v4df vr2_1 =  vr2_1x + vr2_1y + vr2_1z;

    v4df vr2_2x = vdq_a2 * vdq_a2;
    v4df vr2_2y = _mm256_permute4x64_pd(vr2_2x, 201);
    v4df vr2_2z = _mm256_permute4x64_pd(vr2_2x, 210);
    v4df vr2_2 =  vr2_2x + vr2_2y + vr2_2z;

    v4df vr2_3x = vdq_a3 * vdq_a3;
    v4df vr2_3y = _mm256_permute4x64_pd(vr2_3x, 201);
    v4df vr2_3z = _mm256_permute4x64_pd(vr2_3x, 210);
    v4df vr2_3 =  vr2_3x + vr2_3y + vr2_3z;

    v4df vr2_4x = vdq_a4 * vdq_a4;
    v4df vr2_4y = _mm256_permute4x64_pd(vr2_4x, 201);
    v4df vr2_4z = _mm256_permute4x64_pd(vr2_4x, 210);
    v4df vr2_4 =  vr2_4x + vr2_4y + vr2_4z;

    v4df vr2_13 = _mm256_unpacklo_pd(vr2_1, vr2_3);
    v4df vr2_24 = _mm256_unpacklo_pd(vr2_2, vr2_4);
    v4df vr2 = _mm256_shuffle_pd(vr2_13, vr2_24, 12);
    v4df vr6 = vr2 * vr2 * vr2;
    v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);

    v4df mask = vcl2 - vr2;
    vdf = _mm256_blendv_pd(vdf, vzero, mask);

    v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
    v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
    v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
    v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

    v4df vpi_1 = _mm256_load_pd((double*)(p + i_a1));
    vpi_1 += vdq_a1 * vdf_1;
    _mm256_store_pd((double*)(p + i_a1), vpi_1);
    v4df vpj_1 = _mm256_load_pd((double*)(p + j_a1));
    vpj_1 -= vdq_a1 * vdf_1;
    _mm256_store_pd((double*)(p + j_a1), vpj_1);

    v4df vpi_2 = _mm256_load_pd((double*)(p + i_a2));
    vpi_2 += vdq_a2 * vdf_2;
    _mm256_store_pd((double*)(p + i_a2), vpi_2);
    v4df vpj_2 = _mm256_load_pd((double*)(p + j_a2));
    vpj_2 -= vdq_a2 * vdf_2;
    _mm256_store_pd((double*)(p + j_a2), vpj_2);

    v4df vpi_3 = _mm256_load_pd((double*)(p + i_a3));
    vpi_3 += vdq_a3 * vdf_3;
    _mm256_store_pd((double*)(p + i_a3), vpi_3);
    v4df vpj_3 = _mm256_load_pd((double*)(p + j_a3));
    vpj_3 -= vdq_a3 * vdf_3;
    _mm256_store_pd((double*)(p + j_a3), vpj_3);

    v4df vpi_4 = _mm256_load_pd((double*)(p + i_a4));
    vpi_4 += vdq_a4 * vdf_4;
    _mm256_store_pd((double*)(p + i_a4), vpi_4);
    v4df vpj_4 = _mm256_load_pd((double*)(p + j_a4));
    vpj_4 -= vdq_a4 * vdf_4;
    _mm256_store_pd((double*)(p + j_a4), vpj_4);

    i_a1 = i_b1;
    j_a1 = j_b1;
    i_a2 = i_b2;
    j_a2 = j_b2;
    i_a3 = i_b3;
    j_a3 = j_b3;
    i_a4 = i_b4;
    j_a4 = j_b4;
  }
  for (k = (number_of_pairs) / 4 * 4 - 4; k < number_of_pairs; k++) {
    const int i = i_particles[k];
    const int j = j_particles[k];
    double dx = q[j][X] - q[i][X];
    double dy = q[j][Y] - q[i][Y];
    double dz = q[j][Z] - q[i][Z];
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    p[i][X] += df * dx;
    p[i][Y] += df * dy;
    p[i][Z] += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
force_sorted(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      //if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df=0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
void
force_sorted_swp(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    int ja = sorted_list[kp];
    double dxa = q[ja][X] - qx_key;
    double dya = q[ja][Y] - qy_key;
    double dza = q[ja][Z] - qz_key;
    double df = 0.0;
    double dxb = 0.0, dyb = 0.0, dzb = 0.0;
    int jb = 0;

    const int np = number_of_partners[i];
    for (int k = kp; k < np + kp; k++) {
      const double dx = dxa;
      const double dy = dya;
      const double dz = dza;
      double r2 = (dx * dx + dy * dy + dz * dz);
      const int j = ja;
      ja = sorted_list[k + 1];
      dxa = q[ja][X] - qx_key;
      dya = q[ja][Y] - qy_key;
      dza = q[ja][Z] - qz_key;
      if (r2 > CL2)continue;
      pfx += df * dxb;
      pfy += df * dyb;
      pfz += df * dzb;
      p[jb][X] -= df * dxb;
      p[jb][Y] -= df * dyb;
      p[jb][Z] -= df * dzb;
      const double r6 = r2 * r2 * r2;
      df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      jb = j;
      dxb = dx;
      dyb = dy;
      dzb = dz;
    }
    p[jb][X] -= df * dxb;
    p[jb][Y] -= df * dyb;
    p[jb][Z] -= df * dzb;
    p[i][X] += pfx + df * dxb;
    p[i][Y] += pfy + df * dyb;
    p[i][Z] += pfz + df * dzb;
  }
}
//----------------------------------------------------------------------
void
force_sorted_swp_intrin(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    int ja_1 = sorted_list[kp];
    int ja_2 = sorted_list[kp + 1];
    int ja_3 = sorted_list[kp + 2];
    int ja_4 = sorted_list[kp + 3];
    double dxa_1 = q[ja_1][X] - qx_key;
    double dya_1 = q[ja_1][Y] - qy_key;
    double dza_1 = q[ja_1][Z] - qz_key;
    double dxa_2 = q[ja_2][X] - qx_key;
    double dya_2 = q[ja_2][Y] - qy_key;
    double dza_2 = q[ja_2][Z] - qz_key;
    double dxa_3 = q[ja_3][X] - qx_key;
    double dya_3 = q[ja_3][Y] - qy_key;
    double dza_3 = q[ja_3][Z] - qz_key;
    double dxa_4 = q[ja_4][X] - qx_key;
    double dya_4 = q[ja_4][Y] - qy_key;
    double dza_4 = q[ja_4][Z] - qz_key;

    double df_1 = 0.0 , df_2 = 0.0, df_3 = 0.0, df_4 = 0.0;
    double dxb_1 = 0.0, dyb_1 = 0.0, dzb_1 = 0.0;
    double dxb_2 = 0.0, dyb_2 = 0.0, dzb_2 = 0.0;
    double dxb_3 = 0.0, dyb_3 = 0.0, dzb_3 = 0.0;
    double dxb_4 = 0.0, dyb_4 = 0.0, dzb_4 = 0.0;
    int jb_1 = 0, jb_2 = 0, jb_3 = 0, jb_4 = 0;
    const int np = number_of_partners[i];
    for (int k = 0; k < (np/4)*4; k+=4) {
      const double dx_1 = dxa_1;
      const double dy_1 = dya_1;
      const double dz_1 = dza_1;
      double r2_1 = (dx_1 * dx_1 + dy_1 * dy_1 + dz_1 * dz_1);
      const int j_1 = ja_1;
      ja_1 = sorted_list[kp + k + 4];
      dxa_1 = q[ja_1][X] - qx_key;
      dya_1 = q[ja_1][Y] - qy_key;
      dza_1 = q[ja_1][Z] - qz_key;
      pfx += df_1 * dxb_1;
      pfy += df_1 * dyb_1;
      pfz += df_1 * dzb_1;
      p[jb_1][X] -= df_1 * dxb_1;
      p[jb_1][Y] -= df_1 * dyb_1;
      p[jb_1][Z] -= df_1 * dzb_1;
      const double r6_1 = r2_1 * r2_1 * r2_1;
      df_1 = ((24.0 * r6_1 - 48.0) / (r6_1 * r6_1 * r2_1)) * dt;
      if (r2_1 > CL2) df_1 = 0.0;
      jb_1 = j_1;
      dxb_1 = dx_1;
      dyb_1 = dy_1;
      dzb_1 = dz_1;

      const double dx_2 = dxa_2;
      const double dy_2 = dya_2;
      const double dz_2 = dza_2;
      double r2_2 = (dx_2 * dx_2 + dy_2 * dy_2 + dz_2 * dz_2);
      const int j_2 = ja_2;
      ja_2 = sorted_list[kp + k + 5];
      dxa_2 = q[ja_2][X] - qx_key;
      dya_2 = q[ja_2][Y] - qy_key;
      dza_2 = q[ja_2][Z] - qz_key;
      pfx += df_2 * dxb_2;
      pfy += df_2 * dyb_2;
      pfz += df_2 * dzb_2;
      p[jb_2][X] -= df_2 * dxb_2;
      p[jb_2][Y] -= df_2 * dyb_2;
      p[jb_2][Z] -= df_2 * dzb_2;
      const double r6_2 = r2_2 * r2_2 * r2_2;
      df_2 = ((24.0 * r6_2 - 48.0) / (r6_2 * r6_2 * r2_2)) * dt;
      if (r2_2 > CL2) df_2 = 0.0;
      jb_2 = j_2;
      dxb_2 = dx_2;
      dyb_2 = dy_2;
      dzb_2 = dz_2;

      const double dx_3 = dxa_3;
      const double dy_3 = dya_3;
      const double dz_3 = dza_3;
      double r2_3 = (dx_3 * dx_3 + dy_3 * dy_3 + dz_3 * dz_3);
      const int j_3 = ja_3;
      ja_3 = sorted_list[kp + k + 6];
      dxa_3 = q[ja_3][X] - qx_key;
      dya_3 = q[ja_3][Y] - qy_key;
      dza_3 = q[ja_3][Z] - qz_key;
      pfx += df_3 * dxb_3;
      pfy += df_3 * dyb_3;
      pfz += df_3 * dzb_3;
      p[jb_3][X] -= df_3 * dxb_3;
      p[jb_3][Y] -= df_3 * dyb_3;
      p[jb_3][Z] -= df_3 * dzb_3;
      const double r6_3 = r2_3 * r2_3 * r2_3;
      df_3 = ((24.0 * r6_3 - 48.0) / (r6_3 * r6_3 * r2_3)) * dt;
      if (r2_3 > CL2) df_3 = 0.0;
      jb_3 = j_3;
      dxb_3 = dx_3;
      dyb_3 = dy_3;
      dzb_3 = dz_3;

      const double dx_4 = dxa_4;
      const double dy_4 = dya_4;
      const double dz_4 = dza_4;
      double r2_4 = (dx_4 * dx_4 + dy_4 * dy_4 + dz_4 * dz_4);
      const int j_4 = ja_4;
      ja_4 = sorted_list[kp + k + 7];
      dxa_4 = q[ja_4][X] - qx_key;
      dya_4 = q[ja_4][Y] - qy_key;
      dza_4 = q[ja_4][Z] - qz_key;
      pfx += df_4 * dxb_4;
      pfy += df_4 * dyb_4;
      pfz += df_4 * dzb_4;
      p[jb_4][X] -= df_4 * dxb_4;
      p[jb_4][Y] -= df_4 * dyb_4;
      p[jb_4][Z] -= df_4 * dzb_4;
      const double r6_4 = r2_4 * r2_4 * r2_4;
      df_4 = ((24.0 * r6_4 - 48.0) / (r6_4 * r6_4 * r2_4)) * dt;
      if (r2_4 > CL2) df_4 = 0.0;
      jb_4 = j_4;
      dxb_4 = dx_4;
      dyb_4 = dy_4;
      dzb_4 = dz_4;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;

    p[jb_1][X] -= df_1 * dxb_1;
    p[jb_1][Y] -= df_1 * dyb_1;
    p[jb_1][Z] -= df_1 * dzb_1;
    p[i][X] += df_1 * dxb_1;
    p[i][Y] += df_1 * dyb_1;
    p[i][Z] += df_1 * dzb_1;

    p[jb_2][X] -= df_2 * dxb_2;
    p[jb_2][Y] -= df_2 * dyb_2;
    p[jb_2][Z] -= df_2 * dzb_2;
    p[i][X] += df_2 * dxb_2;
    p[i][Y] += df_2 * dyb_2;
    p[i][Z] += df_2 * dzb_2;

    p[jb_3][X] -= df_3 * dxb_3;
    p[jb_3][Y] -= df_3 * dyb_3;
    p[jb_3][Z] -= df_3 * dzb_3;
    p[i][X] += df_3 * dxb_3;
    p[i][Y] += df_3 * dyb_3;
    p[i][Z] += df_3 * dzb_3;

    p[jb_4][X] -= df_4 * dxb_4;
    p[jb_4][Y] -= df_4 * dyb_4;
    p[jb_4][Z] -= df_4 * dzb_4;
    p[i][X] += df_4 * dxb_4;
    p[i][Y] += df_4 * dyb_4;
    p[i][Z] += df_4 * dzb_4;
    pfx = 0.0;
    pfy = 0.0;
    pfz = 0.0;
    for (int k = (np/4)*4; k < np; k++) {
      const int j = sorted_list[k+kp];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df=0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
void
force_sorted_intrin(void) {
  const v4df vzero = _mm256_set_pd(0, 0, 0, 0);
  const v4df vcl2 = _mm256_set_pd(CL2, CL2, CL2, CL2);
  const v4df vc24 = _mm256_set_pd(24 * dt, 24 * dt, 24 * dt, 24 * dt);
  const v4df vc48 = _mm256_set_pd(48 * dt, 48 * dt, 48 * dt, 48 * dt);
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(q + i));
    v4df vpi = _mm256_load_pd((double*)(p + i));
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_a = sorted_list[kp + k];
      v4df vqj_a = _mm256_load_pd((double*)(q + j_a));
      v4df vdq_a = (vqj_a - vqi);
      v4df vd1_a = vdq_a * vdq_a;
      v4df vd2_a = _mm256_permute4x64_pd(vd1_a, 201);
      v4df vd3_a = _mm256_permute4x64_pd(vd1_a, 210);
      v4df vr2_a = vd1_a + vd2_a + vd3_a;

      const int j_b = sorted_list[kp + k + 1];
      v4df vqj_b = _mm256_load_pd((double*)(q + j_b));
      v4df vdq_b = (vqj_b - vqi);
      v4df vd1_b = vdq_b * vdq_b;
      v4df vd2_b = _mm256_permute4x64_pd(vd1_b, 201);
      v4df vd3_b = _mm256_permute4x64_pd(vd1_b, 210);
      v4df vr2_b = vd1_b + vd2_b + vd3_b;

      const int j_c = sorted_list[kp + k + 2];
      v4df vqj_c = _mm256_load_pd((double*)(q + j_c));
      v4df vdq_c = (vqj_c - vqi);
      v4df vd1_c = vdq_c * vdq_c;
      v4df vd2_c = _mm256_permute4x64_pd(vd1_c, 201);
      v4df vd3_c = _mm256_permute4x64_pd(vd1_c, 210);
      v4df vr2_c = vd1_c + vd2_c + vd3_c;

      const int j_d = sorted_list[kp + k + 3];
      v4df vqj_d = _mm256_load_pd((double*)(q + j_d));
      v4df vdq_d = (vqj_d - vqi);
      v4df vd1_d = vdq_d * vdq_d;
      v4df vd2_d = _mm256_permute4x64_pd(vd1_d, 201);
      v4df vd3_d = _mm256_permute4x64_pd(vd1_d, 210);
      v4df vr2_d = vd1_d + vd2_d + vd3_d;

      v4df vr2_ac = _mm256_unpacklo_pd(vr2_a, vr2_c);
      v4df vr2_bd = _mm256_unpacklo_pd(vr2_b, vr2_d);
      v4df vr2 = _mm256_shuffle_pd(vr2_ac, vr2_bd, 12);

      v4df vr6 = vr2 * vr2 * vr2;
      v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      v4df vdf_a = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_b = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_c = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_d = _mm256_permute4x64_pd(vdf, 255);

      v4df vpj_a = _mm256_load_pd((double*)(p + j_a));
      vpi += vdq_a * vdf_a;
      vpj_a -= vdq_a * vdf_a;
      _mm256_store_pd((double*)(p + j_a), vpj_a);

      v4df vpj_b = _mm256_load_pd((double*)(p + j_b));
      vpi += vdq_b * vdf_b;
      vpj_b -= vdq_b * vdf_b;
      _mm256_store_pd((double*)(p + j_b), vpj_b);

      v4df vpj_c = _mm256_load_pd((double*)(p + j_c));
      vpi += vdq_c * vdf_c;
      vpj_c -= vdq_c * vdf_c;
      _mm256_store_pd((double*)(p + j_c), vpj_c);

      v4df vpj_d = _mm256_load_pd((double*)(p + j_d));
      vpi += vdq_d * vdf_d;
      vpj_d -= vdq_d * vdf_d;
      _mm256_store_pd((double*)(p + j_d), vpj_d);
    }
    _mm256_store_pd((double*)(p + i), vpi);
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - q[i][X];
      double dy = q[j][Y] - q[i][Y];
      double dz = q[j][Z] - q[i][Z];
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      p[i][X] += df * dx;
      p[i][Y] += df * dy;
      p[i][Z] += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
  }
}
//----------------------------------------------------------------------
void
measure(void(*pfunc)(), const char *name) {
  double st = myclock();
  const int LOOP = 100;
  for (int i = 0; i < LOOP; i++) {
    pfunc();
  }
  double t = myclock() - st;
  fprintf(stderr, "N=%d, %s %f [sec]\n", particle_number, name, t);
}
//----------------------------------------------------------------------
void
loadpair(void) {
  std::ifstream ifs("pair.dat", std::ios::binary);
  ifs.read((char*)&number_of_pairs, sizeof(int));
  ifs.read((char*)number_of_partners, sizeof(int)*N);
  ifs.read((char*)i_particles, sizeof(int)*MAX_PAIRS);
  ifs.read((char*)j_particles, sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
savepair(void) {
  makepair();
  std::ofstream ofs("pair.dat", std::ios::binary);
  ofs.write((char*)&number_of_pairs, sizeof(int));
  ofs.write((char*)number_of_partners, sizeof(int)*N);
  ofs.write((char*)i_particles, sizeof(int)*MAX_PAIRS);
  ofs.write((char*)j_particles, sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
print_result(void) {
  for (int i = 0; i < 5; i++) {
    printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
  }
  for (int i = particle_number - 5; i < particle_number; i++) {
    printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
  }
}
//----------------------------------------------------------------------
int
main(void) {
  init();
  struct stat st;
  int ret = stat("pair.dat", &st);
  if (ret == 0) {
    std::cerr << "A pair-file is found. I use it." << std::endl;
    loadpair();
  } else {
    std::cerr << "Make pairlist." << std::endl;
    savepair();
  }
  std::cerr << "Number of pairs: " << number_of_pairs << std::endl;
  sortpair();
#ifdef PAIR
  measure(&force_pair, "pair");
  print_result();
#elif SORTED
  measure(&force_sorted, "sorted");
  print_result();
#elif S_INTRIN
  measure(&force_sorted_intrin, "sorted_intrin");
  print_result();
#elif S_SWP
  measure(&force_sorted_swp, "sorted_swp");
  print_result();
#elif S_SWP_INTRIN
  measure(&force_sorted_swp_intrin, "sorted_swp_intrin");
  print_result();
#elif P_SWP
  measure(&force_pair_swp, "pair_swp");
  print_result();
#elif P_INTRIN
  measure(&force_pair_intrin, "pair_intrin");
  print_result();
#else
  measure(&force_pair, "pair");
  measure(&force_pair_swp, "pair_swp");
  measure(&force_pair_intrin, "pair_intrin");
  measure(&force_sorted, "sorted");
  measure(&force_sorted_swp, "sorted_swp");
  measure(&force_sorted_intrin, "sorted_intrin");
  measure(&force_sorted_swp_intrin, "sorted_swp_intrin");
#endif
}
//----------------------------------------------------------------------

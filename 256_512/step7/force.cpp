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
const int N = 400000;
const int MAX_PAIRS = 30 * N;
double L = 50.0;
const double dt = 0.001;
const int D = 4;
enum {X = 0, Y, Z, W, PX, PY, PZ, WX};
double q[N][D];
double p[N][D];
__attribute__((aligned(64))) double z[N][8];

double q2[D][N];
double p2[D][N];

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
typedef double v8df __attribute__((vector_size(64)));
//----------------------------------------------------------------------
void
copy_to_z(void){
  for(int i=0;i<particle_number;i++){
    z[i][X] = q[i][X];
    z[i][Y] = q[i][Y];
    z[i][Z] = q[i][Z];
    z[i][PX] = p[i][X];
    z[i][PY] = p[i][Y];
    z[i][PZ] = p[i][Z];
  }
}
//----------------------------------------------------------------------
void
copy_to_2(void){
  for(int i=0;i<particle_number;i++){
    q2[X][i] = q[i][X];
    q2[Y][i] = q[i][Y];
    q2[Z][i] = q[i][Z];
    p2[X][i] = p[i][X];
    p2[Y][i] = p[i][Y];
    p2[Z][i] = p[i][Z];
  }
}
//----------------------------------------------------------------------
void
copy_from_2(void){
  for(int i=0;i<particle_number;i++){
    q[i][X] = q2[X][i];
    q[i][Y] = q2[Y][i];
    q[i][Z] = q2[Z][i];
    p[i][X] = p2[X][i];
    p[i][Y] = p2[Y][i];
    p[i][Z] = p2[Z][i];
  }
}

//----------------------------------------------------------------------
void
copy_from_z(void){
  for(int i=0;i<particle_number;i++){
    q[i][X] = z[i][X];
    q[i][Y] = z[i][Y];
    q[i][Z] = z[i][Z];
    p[i][X] = z[i][PX];
    p[i][Y] = z[i][PY];
    p[i][Z] = z[i][PZ];
  }
}
//----------------------------------------------------------------------
template <class T>
void puts(T r) {
  double *a = (double*)(&r);
  for(int i=0;i < sizeof(T)/sizeof(double);i++){
    printf("%.10f ",a[i]);
  }
  printf("\n");
}
//----------------------------------------------------------------------
void puts4(v4df r) {
  double *a = (double*)(&r);
  for(int i=0;i < sizeof(v4df)/sizeof(double);i++){
    printf("%.10f ",a[i]);
  }
  printf("\n");
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
  std::mt19937 mt(5);
  std::uniform_real_distribution<double> ud(0.0,1.0);
  for (int i = 0; i < particle_number; i++) {
    p[i][X] = 0.0;
    p[i][Y] = 0.0;
    p[i][Z] = 0.0;
    //p[i][X] = ud(mt);
    //p[i][Y] = ud(mt);
    //p[i][Z] = ud(mt);
  }
}
//----------------------------------------------------------------------
__attribute__((noinline))
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
__attribute__((noinline))
void
force_pair_z(void) {
  for (int k = 0; k < number_of_pairs; k++) {
    const int i = i_particles[k];
    const int j = j_particles[k];
    double dx = z[j][X] - z[i][X];
    double dy = z[j][Y] - z[i][Y];
    double dz = z[j][Z] - z[i][Z];
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    z[i][PX] += df * dx;
    z[i][PY] += df * dy;
    z[i][PZ] += df * dz;
    z[j][PX] -= df * dx;
    z[j][PY] -= df * dy;
    z[j][PZ] -= df * dz;
  }
}
//----------------------------------------------------------------------
__attribute__((noinline))
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
#define p4(x) printf("%.10f %.10f %.10f %.10f\n",x##_1,x##_2,x##_3,x##_4);
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_pair_swp_intrin(void) {
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
__attribute__((noinline))
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
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
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
__attribute__((noinline))
void
force_sorted2(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q2[X][i];
    const double qy_key = q2[Y][i];
    const double qz_key = q2[Z][i];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q2[X][j] - qx_key;
      double dy = q2[Y][j] - qy_key;
      double dz = q2[Z][j] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p2[X][j] -= df * dx;
      p2[Y][j] -= df * dy;
      p2[Z][j] -= df * dz;
    }
    p2[X][i] += pfx;
    p2[Y][i] += pfy;
    p2[Z][i] += pfz;
  }
}
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_sorted_z(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = z[i][X];
    const double qy_key = z[i][Y];
    const double qz_key = z[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = z[j][X] - qx_key;
      double dy = z[j][Y] - qy_key;
      double dz = z[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df=0.0; 
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      z[j][PX] -= df * dx;
      z[j][PY] -= df * dy;
      z[j][PZ] -= df * dz;
    }
    z[i][PX] += pfx;
    z[i][PY] += pfy;
    z[i][PZ] += pfz;
  }
}
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_sorted_z_intrin(void) {
  const int pn = particle_number;
  const v8df vzero = _mm512_setzero_pd();
  const v8df vcl2 = _mm512_set1_pd(CL2);
  const v8df vc24 = _mm512_set1_pd(24.0*dt);
  const v8df vc48 = _mm512_set1_pd(48.0*dt);
  const __m512i idx = _mm512_set_epi64(3,2,1,0,7,6,5,4);
  const __m512i idx_q = _mm512_set_epi64(11,10,9,8,3,2,1,0);
  const __m512i idx_p = _mm512_set_epi64(15,14,13,12,7,6,5,4);
  const __m512i idx_f12 = _mm512_set_epi64(1,1,1,1,0,0,0,0);
  const __m512i idx_f34 = _mm512_set_epi64(3,3,3,3,2,2,2,2);
  const __m512i idx_f56 = _mm512_set_epi64(5,5,5,5,4,4,4,4);
  const __m512i idx_f78 = _mm512_set_epi64(7,7,7,7,6,6,6,6);
  const __m512i idx_0123 = _mm512_set_epi64(3,2,1,0,3,2,1,0);
  const __m512i idx2 = _mm512_set_epi64(13,9,12,8,5,1,4,0);
  const __m512i idx3 = _mm512_set_epi64(15,11,14,10,7,3,6,2);
  const __mmask8 khigh = 16+32+64+128;
  for (int i = 0; i < pn; i++) {
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    v8df vzi = _mm512_loadu_pd((double*)(z+i));
    v8df vqi = _mm512_permutexvar_pd(idx_0123, vzi);;
    v8df vpi = _mm512_setzero_pd();
    for (int k = 0; k < (np/8*8); k+=8) {
      const int j_1 = sorted_list[kp + k]; 
      const int j_2 = sorted_list[kp + k + 1]; 
      const int j_3 = sorted_list[kp + k + 2]; 
      const int j_4 = sorted_list[kp + k + 3]; 
      const int j_5 = sorted_list[kp + k + 4]; 
      const int j_6 = sorted_list[kp + k + 5]; 
      const int j_7 = sorted_list[kp + k + 6]; 
      const int j_8 = sorted_list[kp + k + 7]; 
      v8df vzj_1 = _mm512_loadu_pd((double*)(z+j_1));
      v8df vzj_2 = _mm512_loadu_pd((double*)(z+j_2));
      v8df vqj_12= _mm512_permutex2var_pd(vzj_1, idx_q, vzj_2);
      v8df vpj_12= _mm512_permutex2var_pd(vzj_1, idx_p, vzj_2);
      v8df vdq_12 = vqj_12 - vqi;

      v8df vzj_3 = _mm512_loadu_pd((double*)(z+j_3));
      v8df vzj_4 = _mm512_loadu_pd((double*)(z+j_4));
      v8df vqj_34= _mm512_permutex2var_pd(vzj_3, idx_q, vzj_4);
      v8df vpj_34= _mm512_permutex2var_pd(vzj_3, idx_p, vzj_4);
      v8df vdq_34 = vqj_34 - vqi;

      v8df vzj_5 = _mm512_loadu_pd((double*)(z+j_5));
      v8df vzj_6 = _mm512_loadu_pd((double*)(z+j_6));
      v8df vqj_56= _mm512_permutex2var_pd(vzj_5, idx_q, vzj_6);
      v8df vpj_56= _mm512_permutex2var_pd(vzj_5, idx_p, vzj_6);
      v8df vdq_56 = vqj_56 - vqi;

      v8df vzj_7 = _mm512_loadu_pd((double*)(z+j_7));
      v8df vzj_8 = _mm512_loadu_pd((double*)(z+j_8));
      v8df vqj_78= _mm512_permutex2var_pd(vzj_7, idx_q, vzj_8);
      v8df vpj_78= _mm512_permutex2var_pd(vzj_7, idx_p, vzj_8);
      v8df vdq_78 = vqj_78 - vqi;

      v8df tmp0 = _mm512_unpacklo_pd(vdq_12, vdq_34);
      v8df tmp1 = _mm512_unpackhi_pd(vdq_12, vdq_34);
      v8df tmp2 = _mm512_unpacklo_pd(vdq_56, vdq_78);
      v8df tmp3 = _mm512_unpackhi_pd(vdq_56, vdq_78);

      v8df vdx = _mm512_permutex2var_pd(tmp0, idx2, tmp2);
      v8df vdy = _mm512_permutex2var_pd(tmp1, idx2, tmp3);
      v8df vdz = _mm512_permutex2var_pd(tmp0, idx3, tmp2);

      v8df vr2 = vdx*vdx + vdy*vdy + vdz*vdz;
      v8df vr6 = vr2 * vr2 * vr2;
      v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);

      __mmask8 kcmp = _mm512_cmp_pd_mask(vr2,vcl2, _CMP_GT_OS);
      vdf = _mm512_mask_blend_pd(kcmp, vdf, vzero);
      v8df vdf_12 = _mm512_permutexvar_pd(idx_f12, vdf);
      v8df vdf_34 = _mm512_permutexvar_pd(idx_f34, vdf);
      v8df vdf_56 = _mm512_permutexvar_pd(idx_f56, vdf);
      v8df vdf_78 = _mm512_permutexvar_pd(idx_f78, vdf);
      vpj_12 -= vdf_12 * vdq_12;
      vpj_34 -= vdf_34 * vdq_34;
      vpj_56 -= vdf_56 * vdq_56;
      vpj_78 -= vdf_78 * vdq_78;
      vpi += vdf_12 * vdq_12;
      vpi += vdf_34 * vdq_34;
      vpi += vdf_56 * vdq_56;
      vpi += vdf_78 * vdq_78;
      v8df vpj_11 = _mm512_permutexvar_pd(idx_0123, vpj_12);
      v8df vpj_33 = _mm512_permutexvar_pd(idx_0123, vpj_34);
      v8df vpj_55 = _mm512_permutexvar_pd(idx_0123, vpj_56);
      v8df vpj_77 = _mm512_permutexvar_pd(idx_0123, vpj_78);
      _mm512_mask_store_pd((double*)(z+j_1), khigh, vpj_11);
      _mm512_mask_store_pd((double*)(z+j_3), khigh, vpj_33);
      _mm512_mask_store_pd((double*)(z+j_5), khigh, vpj_55);
      _mm512_mask_store_pd((double*)(z+j_7), khigh, vpj_77);
      _mm512_mask_store_pd((double*)(z+j_2), khigh, vpj_12);
      _mm512_mask_store_pd((double*)(z+j_4), khigh, vpj_34);
      _mm512_mask_store_pd((double*)(z+j_6), khigh, vpj_56);
      _mm512_mask_store_pd((double*)(z+j_8), khigh, vpj_78);
    }
    v4df vpi_low = _mm512_extractf64x4_pd(vpi, 0);
    v4df vpi_high = _mm512_extractf64x4_pd(vpi, 1);
    v4df vdpi = vpi_low + vpi_high;
    v8df vzdpi = _mm512_insertf64x4(vzero, vdpi, 1); 
    vzi += vzdpi;
    _mm512_store_pd((double*)(z+i), vzi);
    const double qix = z[i][X];
    const double qiy = z[i][Y];
    const double qiz = z[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    for (int k = (np/8*8); k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = z[j][X] - qix;
      double dy = z[j][Y] - qiy;
      double dz = z[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df=0.0; 
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      z[j][PX] -= df * dx;
      z[j][PY] -= df * dy;
      z[j][PZ] -= df * dz;
    }
    z[i][PX] += pfx;
    z[i][PY] += pfy;
    z[i][PZ] += pfz;
  }
}
//----------------------------------------------------------------------
// NOTE: gather/scatter + swp
__attribute__((noinline))
void
force_sorted_z_intrin_gs_swp(void) {
  const auto vc24  = _mm512_set1_pd(24.0 * dt);
  const auto vc48  = _mm512_set1_pd(48.0 * dt);
  const auto vcl2  = _mm512_set1_pd(CL2);
  const auto v2    = _mm512_set1_pd(2.0);
  const auto vzero = _mm512_setzero_pd();
  const auto pn = particle_number;
  const auto vpitch = _mm512_set1_epi64(8);

  for (int i = 0; i < pn; i++) {
    const auto vqxi = _mm512_set1_pd(z[i][X]);
    const auto vqyi = _mm512_set1_pd(z[i][Y]);
    const auto vqzi = _mm512_set1_pd(z[i][Z]);

    auto vpxi = _mm512_setzero_pd();
    auto vpyi = _mm512_setzero_pd();
    auto vpzi = _mm512_setzero_pd();

    const auto np = number_of_partners[i];
    const auto kp = pointer[i];
    const int* ptr_list = &sorted_list[kp];

    const auto vnp = _mm512_set1_epi64(np);
    auto vk_idx = _mm512_set_epi64(7LL, 6LL, 5LL, 4LL,
                                   3LL, 2LL, 1LL, 0LL);

    // initial force calculation
    // load position
    auto vindex_a = _mm256_slli_epi32(_mm256_lddqu_si256((const __m256i*)ptr_list), 3);
    auto mask_a = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
    auto vqxj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][X], 8);
    auto vqyj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][Y], 8);
    auto vqzj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][Z], 8);

    // calc distance
    auto vdx_a = _mm512_sub_pd(vqxj, vqxi);
    auto vdy_a = _mm512_sub_pd(vqyj, vqyi);
    auto vdz_a = _mm512_sub_pd(vqzj, vqzi);
    auto vr2 = _mm512_fmadd_pd(vdz_a,
                               vdz_a,
                               _mm512_fmadd_pd(vdy_a,
                                               vdy_a,
                                               _mm512_mul_pd(vdx_a, vdx_a)));

    // calc force norm
    auto vr6 = _mm512_mul_pd(_mm512_mul_pd(vr2, vr2), vr2);
    auto vdf = _mm512_div_pd(_mm512_fmsub_pd(vc24, vr6, vc48),
                             _mm512_mul_pd(_mm512_mul_pd(vr6, vr6), vr2));

    auto le_cl2 = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
    mask_a = _mm512_kand(mask_a, le_cl2);
    vdf = _mm512_mask_blend_pd(mask_a, vzero, vdf);

    for (int k = 8; k < np; k += 8) {
      // load position
      ptr_list += 8;
      auto vindex_b = _mm256_slli_epi32(_mm256_lddqu_si256((const __m256i*)ptr_list), 3);
      vk_idx = _mm512_add_epi64(vk_idx, vpitch);
      auto mask_b = _mm512_cmp_epi64_mask(vk_idx,
                                          vnp,
                                          _MM_CMPINT_LT);
      vqxj = _mm512_mask_i32gather_pd(vzero, mask_b, vindex_b, &z[0][X], 8);
      vqyj = _mm512_mask_i32gather_pd(vzero, mask_b, vindex_b, &z[0][Y], 8);
      vqzj = _mm512_mask_i32gather_pd(vzero, mask_b, vindex_b, &z[0][Z], 8);

      // calc distance
      auto vdx_b = _mm512_sub_pd(vqxj, vqxi);
      auto vdy_b = _mm512_sub_pd(vqyj, vqyi);
      auto vdz_b = _mm512_sub_pd(vqzj, vqzi);
      vr2 = _mm512_fmadd_pd(vdz_b,
                            vdz_b,
                            _mm512_fmadd_pd(vdy_b,
                                            vdy_b,
                                            _mm512_mul_pd(vdx_b,
                                                          vdx_b)));

      // write back j particle momentum
      vpxi = _mm512_fmadd_pd(vdf, vdx_a, vpxi);
      vpyi = _mm512_fmadd_pd(vdf, vdy_a, vpyi);
      vpzi = _mm512_fmadd_pd(vdf, vdz_a, vpzi);

      auto vpxj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][PX], 8);
      auto vpyj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][PY], 8);
      auto vpzj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][PZ], 8);

      vpxj = _mm512_fnmadd_pd(vdf, vdx_a, vpxj);
      vpyj = _mm512_fnmadd_pd(vdf, vdy_a, vpyj);
      vpzj = _mm512_fnmadd_pd(vdf, vdz_a, vpzj);

      _mm512_mask_i32scatter_pd(&z[0][PX], mask_a, vindex_a, vpxj, 8);
      _mm512_mask_i32scatter_pd(&z[0][PY], mask_a, vindex_a, vpyj, 8);
      _mm512_mask_i32scatter_pd(&z[0][PZ], mask_a, vindex_a, vpzj, 8);

      // calc force norm
      vr6 = _mm512_mul_pd(_mm512_mul_pd(vr2, vr2), vr2);
      vdf = _mm512_div_pd(_mm512_fmsub_pd(vc24, vr6, vc48),
                          _mm512_mul_pd(_mm512_mul_pd(vr6, vr6), vr2));

      le_cl2 = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      mask_b = _mm512_kand(mask_b, le_cl2);
      vdf = _mm512_mask_blend_pd(mask_b, vzero, vdf);

      // send to next
      vindex_a = vindex_b;
      mask_a   = mask_b;
      vdx_a    = vdx_b;
      vdy_a    = vdy_b;
      vdz_a    = vdz_b;
    } // end of k loop

    // final write back momentum
    // write back j particle momentum
    vpxi = _mm512_fmadd_pd(vdf, vdx_a, vpxi);
    vpyi = _mm512_fmadd_pd(vdf, vdy_a, vpyi);
    vpzi = _mm512_fmadd_pd(vdf, vdz_a, vpzi);

    auto vpxj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][PX], 8);
    auto vpyj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][PY], 8);
    auto vpzj = _mm512_mask_i32gather_pd(vzero, mask_a, vindex_a, &z[0][PZ], 8);

    vpxj = _mm512_fnmadd_pd(vdf, vdx_a, vpxj);
    vpyj = _mm512_fnmadd_pd(vdf, vdy_a, vpyj);
    vpzj = _mm512_fnmadd_pd(vdf, vdz_a, vpzj);

    _mm512_mask_i32scatter_pd(&z[0][PX], mask_a, vindex_a, vpxj, 8);
    _mm512_mask_i32scatter_pd(&z[0][PY], mask_a, vindex_a, vpyj, 8);
    _mm512_mask_i32scatter_pd(&z[0][PZ], mask_a, vindex_a, vpzj, 8);

    // write back i particle momentum
    z[i][PX] += _mm512_reduce_add_pd(vpxi);
    z[i][PY] += _mm512_reduce_add_pd(vpyi);
    z[i][PZ] += _mm512_reduce_add_pd(vpzi);
  } // end of i loop
}
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_sorted_z_intrin_swp(void) {
  const int pn = particle_number;
  const v8df vzero = _mm512_setzero_pd();
  const v8df vcl2 = _mm512_set1_pd(CL2);
  const v8df vc24 = _mm512_set1_pd(24.0*dt);
  const v8df vc48 = _mm512_set1_pd(48.0*dt);
  const __m512i idx = _mm512_set_epi64(3,2,1,0,7,6,5,4);
  const __m512i idx_q = _mm512_set_epi64(11,10,9,8,3,2,1,0);
  const __m512i idx_p = _mm512_set_epi64(15,14,13,12,7,6,5,4);
  const __m512i idx_f12 = _mm512_set_epi64(1,1,1,1,0,0,0,0);
  const __m512i idx_f34 = _mm512_set_epi64(3,3,3,3,2,2,2,2);
  const __m512i idx_f56 = _mm512_set_epi64(5,5,5,5,4,4,4,4);
  const __m512i idx_f78 = _mm512_set_epi64(7,7,7,7,6,6,6,6);
  const __m512i idx_0123 = _mm512_set_epi64(3,2,1,0,3,2,1,0);
  const __m512i idx2 = _mm512_set_epi64(13,9,12,8,5,1,4,0);
  const __m512i idx3 = _mm512_set_epi64(15,11,14,10,7,3,6,2);
  const __mmask8 khigh = 16+32+64+128;
  for (int i = 0; i < pn; i++) {
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    v8df vzi = _mm512_loadu_pd((double*)(z+i));
    v8df vqi = _mm512_permutexvar_pd(idx_0123, vzi);;
    v8df vpi = _mm512_setzero_pd();
    int ja_1 = sorted_list[kp]; 
    int ja_2 = sorted_list[kp + 1]; 
    int ja_3 = sorted_list[kp + 2]; 
    int ja_4 = sorted_list[kp + 3]; 
    int ja_5 = sorted_list[kp + 4]; 
    int ja_6 = sorted_list[kp + 5]; 
    int ja_7 = sorted_list[kp + 6]; 
    int ja_8 = sorted_list[kp + 7]; 
    v8df vzja_1 = _mm512_loadu_pd((double*)(z+ja_1));
    v8df vzja_2 = _mm512_loadu_pd((double*)(z+ja_2));
    v8df vzja_3 = _mm512_loadu_pd((double*)(z+ja_3));
    v8df vzja_4 = _mm512_loadu_pd((double*)(z+ja_4));
    v8df vzja_5 = _mm512_loadu_pd((double*)(z+ja_5));
    v8df vzja_6 = _mm512_loadu_pd((double*)(z+ja_6));
    v8df vzja_7 = _mm512_loadu_pd((double*)(z+ja_7));
    v8df vzja_8 = _mm512_loadu_pd((double*)(z+ja_8));

    v8df vqj_12= _mm512_permutex2var_pd(vzja_1, idx_q, vzja_2);
    v8df vpja_12= _mm512_permutex2var_pd(vzja_1, idx_p, vzja_2);
    v8df vdqa_12 = vqj_12 - vqi;
    v8df vqj_34= _mm512_permutex2var_pd(vzja_3, idx_q, vzja_4);
    v8df vpja_34= _mm512_permutex2var_pd(vzja_3, idx_p, vzja_4);
    v8df vdqa_34 = vqj_34 - vqi;
    v8df vqj_56= _mm512_permutex2var_pd(vzja_5, idx_q, vzja_6);
    v8df vpja_56= _mm512_permutex2var_pd(vzja_5, idx_p, vzja_6);
    v8df vdqa_56 = vqj_56 - vqi;
    v8df vqj_78= _mm512_permutex2var_pd(vzja_7, idx_q, vzja_8);
    v8df vpja_78= _mm512_permutex2var_pd(vzja_7, idx_p, vzja_8);
    v8df vdqa_78 = vqj_78 - vqi;


    for (int k = 0; k < (np/8*8); k+=8) {
      const int j_1 = ja_1;
      const int j_2 = ja_2;
      const int j_3 = ja_3;
      const int j_4 = ja_4;
      const int j_5 = ja_5;
      const int j_6 = ja_6;
      const int j_7 = ja_7;
      const int j_8 = ja_8;
      v8df vzj_1 = vzja_1;
      v8df vzj_2 = vzja_2;
      v8df vzj_3 = vzja_3;
      v8df vzj_4 = vzja_4;
      v8df vzj_5 = vzja_5;
      v8df vzj_6 = vzja_6;
      v8df vzj_7 = vzja_7;
      v8df vzj_8 = vzja_8;
      v8df vpj_12 = vpja_12;
      v8df vpj_34 = vpja_34;
      v8df vpj_56 = vpja_56;
      v8df vpj_78 = vpja_78;
      v8df vdq_12 = vdqa_12;
      v8df vdq_34 = vdqa_34;
      v8df vdq_56 = vdqa_56;
      v8df vdq_78 = vdqa_78;

      ja_1 = sorted_list[kp + k + 8]; 
      ja_2 = sorted_list[kp + k + 1 + 8]; 
      ja_3 = sorted_list[kp + k + 2 + 8]; 
      ja_4 = sorted_list[kp + k + 3 + 8]; 
      ja_5 = sorted_list[kp + k + 4 + 8]; 
      ja_6 = sorted_list[kp + k + 5 + 8]; 
      ja_7 = sorted_list[kp + k + 6 + 8]; 
      ja_8 = sorted_list[kp + k + 7 + 8]; 


      v8df tmp0 = _mm512_unpacklo_pd(vdq_12, vdq_34);
      v8df tmp1 = _mm512_unpackhi_pd(vdq_12, vdq_34);
      v8df tmp2 = _mm512_unpacklo_pd(vdq_56, vdq_78);
      v8df tmp3 = _mm512_unpackhi_pd(vdq_56, vdq_78);

      v8df vdx = _mm512_permutex2var_pd(tmp0, idx2, tmp2);
      v8df vdy = _mm512_permutex2var_pd(tmp1, idx2, tmp3);
      v8df vdz = _mm512_permutex2var_pd(tmp0, idx3, tmp2);

      v8df vr2 = vdx*vdx + vdy*vdy + vdz*vdz;
      v8df vr6 = vr2 * vr2 * vr2;
      v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      vzja_1 = _mm512_loadu_pd((double*)(z+ja_1));
      vzja_2 = _mm512_loadu_pd((double*)(z+ja_2));
      vzja_3 = _mm512_loadu_pd((double*)(z+ja_3));
      vzja_4 = _mm512_loadu_pd((double*)(z+ja_4));
      vzja_5 = _mm512_loadu_pd((double*)(z+ja_5));
      vzja_6 = _mm512_loadu_pd((double*)(z+ja_6));
      vzja_7 = _mm512_loadu_pd((double*)(z+ja_7));
      vzja_8 = _mm512_loadu_pd((double*)(z+ja_8));

      v8df vqja_12= _mm512_permutex2var_pd(vzja_1, idx_q, vzja_2);
      vpja_12= _mm512_permutex2var_pd(vzja_1, idx_p, vzja_2);
      vdqa_12 = vqja_12 - vqi;
      v8df vqj_34= _mm512_permutex2var_pd(vzja_3, idx_q, vzja_4);
      vpja_34= _mm512_permutex2var_pd(vzja_3, idx_p, vzja_4);
      vdqa_34 = vqj_34 - vqi;
      v8df vqj_56= _mm512_permutex2var_pd(vzja_5, idx_q, vzja_6);
      vpja_56= _mm512_permutex2var_pd(vzja_5, idx_p, vzja_6);
      vdqa_56 = vqj_56 - vqi;
      v8df vqj_78= _mm512_permutex2var_pd(vzja_7, idx_q, vzja_8);
      vpja_78= _mm512_permutex2var_pd(vzja_7, idx_p, vzja_8);
      vdqa_78 = vqj_78 - vqi;

      __mmask8 kcmp = _mm512_cmp_pd_mask(vr2,vcl2, _CMP_GT_OS);
      vdf = _mm512_mask_blend_pd(kcmp, vdf, vzero);
      v8df vdf_12 = _mm512_permutexvar_pd(idx_f12, vdf);
      v8df vdf_34 = _mm512_permutexvar_pd(idx_f34, vdf);
      v8df vdf_56 = _mm512_permutexvar_pd(idx_f56, vdf);
      v8df vdf_78 = _mm512_permutexvar_pd(idx_f78, vdf);
      vpj_12 -= vdf_12 * vdq_12;
      vpj_34 -= vdf_34 * vdq_34;
      vpj_56 -= vdf_56 * vdq_56;
      vpj_78 -= vdf_78 * vdq_78;
      vpi += vdf_12 * vdq_12;
      vpi += vdf_34 * vdq_34;
      vpi += vdf_56 * vdq_56;
      vpi += vdf_78 * vdq_78;
      v8df vpj_11 = _mm512_permutexvar_pd(idx_0123, vpj_12);
      v8df vpj_33 = _mm512_permutexvar_pd(idx_0123, vpj_34);
      v8df vpj_55 = _mm512_permutexvar_pd(idx_0123, vpj_56);
      v8df vpj_77 = _mm512_permutexvar_pd(idx_0123, vpj_78);
      _mm512_mask_store_pd((double*)(z+j_1), khigh, vpj_11);
      _mm512_mask_store_pd((double*)(z+j_3), khigh, vpj_33);
      _mm512_mask_store_pd((double*)(z+j_5), khigh, vpj_55);
      _mm512_mask_store_pd((double*)(z+j_7), khigh, vpj_77);
      _mm512_mask_store_pd((double*)(z+j_2), khigh, vpj_12);
      _mm512_mask_store_pd((double*)(z+j_4), khigh, vpj_34);
      _mm512_mask_store_pd((double*)(z+j_6), khigh, vpj_56);
      _mm512_mask_store_pd((double*)(z+j_8), khigh, vpj_78);
    }
    v4df vpi_low = _mm512_extractf64x4_pd(vpi, 0);
    v4df vpi_high = _mm512_extractf64x4_pd(vpi, 1);
    v4df vdpi = vpi_low + vpi_high;
    v8df vzdpi = _mm512_insertf64x4(vzero, vdpi, 1); 
    vzi += vzdpi;
    _mm512_store_pd((double*)(z+i), vzi);
    const double qix = z[i][X];
    const double qiy = z[i][Y];
    const double qiz = z[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    for (int k = (np/8*8); k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = z[j][X] - qix;
      double dy = z[j][Y] - qiy;
      double dz = z[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df=0.0; 
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      z[j][PX] -= df * dx;
      z[j][PY] -= df * dy;
      z[j][PZ] -= df * dz;
    }
    z[i][PX] += pfx;
    z[i][PY] += pfy;
    z[i][PZ] += pfz;
  }

}
//----------------------------------------------------------------------
__attribute__((noinline))
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
#define pv(a,b) printf("%.10f %.10f %.10f %.10f\n",a##x##b,a##y##b,a##z##b,0.0);
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_sorted_swp_intrin(void) {
  const int pn = particle_number;
  const v4df vzero = _mm256_set_pd(0, 0, 0, 0);
  const v4df vcl2 = _mm256_set_pd(CL2, CL2, CL2, CL2);
  const v4df vc24 = _mm256_set_pd(24 * dt, 24 * dt, 24 * dt, 24 * dt);
  const v4df vc48 = _mm256_set_pd(48 * dt, 48 * dt, 48 * dt, 48 * dt);
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(q + i));
    v4df vpf = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    const int kp = pointer[i];
    int ja_1 = sorted_list[kp];
    int ja_2 = sorted_list[kp + 1];
    int ja_3 = sorted_list[kp + 2];
    int ja_4 = sorted_list[kp + 3];
    v4df vqj_1 = _mm256_load_pd((double*)(q + ja_1));
    v4df vdqa_1 = vqj_1 - vqi;
    v4df vqj_2 = _mm256_load_pd((double*)(q + ja_2));
    v4df vdqa_2 = vqj_2 - vqi;
    v4df vqj_3 = _mm256_load_pd((double*)(q + ja_3));
    v4df vdqa_3 = vqj_3 - vqi;
    v4df vqj_4 = _mm256_load_pd((double*)(q + ja_4));
    v4df vdqa_4 = vqj_4 - vqi;

    v4df vdf = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);

    v4df vdqb_1 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_2 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_3 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_4 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);

    int jb_1 = 0, jb_2 = 0, jb_3 = 0, jb_4 = 0;
    const int np = number_of_partners[i];
    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_1 = ja_1;
      const int j_2 = ja_2;
      const int j_3 = ja_3;
      const int j_4 = ja_4;
      v4df vdq_1 = vdqa_1;
      v4df vdq_2 = vdqa_2;
      v4df vdq_3 = vdqa_3;
      v4df vdq_4 = vdqa_4;

      ja_1 = sorted_list[kp + k + 4];
      ja_2 = sorted_list[kp + k + 5];
      ja_3 = sorted_list[kp + k + 6];
      ja_4 = sorted_list[kp + k + 7];

      v4df vr2s_1 = vdq_1 * vdq_1;
      v4df vr2t_1 = _mm256_permute4x64_pd(vr2s_1, 201);
      v4df vr2u_1 = _mm256_permute4x64_pd(vr2s_1, 210);
      v4df vr2_1 = vr2s_1 + vr2t_1 + vr2u_1;

      v4df vr2s_2 = vdq_2 * vdq_2;
      v4df vr2t_2 = _mm256_permute4x64_pd(vr2s_2, 201);
      v4df vr2u_2 = _mm256_permute4x64_pd(vr2s_2, 210);
      v4df vr2_2 = vr2s_2 + vr2t_2 + vr2u_2;

      v4df vr2s_3 = vdq_3 * vdq_3;
      v4df vr2t_3 = _mm256_permute4x64_pd(vr2s_3, 201);
      v4df vr2u_3 = _mm256_permute4x64_pd(vr2s_3, 210);
      v4df vr2_3 = vr2s_3 + vr2t_3 + vr2u_3;

      v4df vr2s_4 = vdq_4 * vdq_4;
      v4df vr2t_4 = _mm256_permute4x64_pd(vr2s_4, 201);
      v4df vr2u_4 = _mm256_permute4x64_pd(vr2s_4, 210);
      v4df vr2_4 = vr2s_4 + vr2t_4 + vr2u_4;

      v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

      vqj_1 = _mm256_load_pd((double*)(q + ja_1));
      vdqa_1 = vqj_1 - vqi;
      vpf += vdf_1 * vdqb_1;

      v4df vpjb_1 = _mm256_load_pd((double*)(p + jb_1));
      vpjb_1 -= vdf_1 * vdqb_1;
      _mm256_store_pd((double*)(p + jb_1), vpjb_1);

      vqj_2 = _mm256_load_pd((double*)(q + ja_2));
      vdqa_2 = vqj_2 - vqi;
      vpf += vdf_2 * vdqb_2;

      v4df vpjb_2 = _mm256_load_pd((double*)(p + jb_2));
      vpjb_2 -= vdf_2 * vdqb_2;
      _mm256_store_pd((double*)(p + jb_2), vpjb_2);

      vqj_3 = _mm256_load_pd((double*)(q + ja_3));
      vdqa_3 = vqj_3 - vqi;
      vpf += vdf_3 * vdqb_3;

      v4df vpjb_3 = _mm256_load_pd((double*)(p + jb_3));
      vpjb_3 -= vdf_3 * vdqb_3;
      _mm256_store_pd((double*)(p + jb_3), vpjb_3);

      vqj_4 = _mm256_load_pd((double*)(q + ja_4));
      vdqa_4 = vqj_4 - vqi;
      vpf += vdf_4 * vdqb_4;

      v4df vpjb_4 = _mm256_load_pd((double*)(p + jb_4));
      vpjb_4 -= vdf_4 * vdqb_4;
      _mm256_store_pd((double*)(p + jb_4), vpjb_4);

      v4df vr2_13 = _mm256_unpacklo_pd(vr2_1, vr2_3);
      v4df vr2_24 = _mm256_unpacklo_pd(vr2_2, vr2_4);
      v4df vr2 = _mm256_shuffle_pd(vr2_13, vr2_24, 12);
      v4df vr6 = vr2 * vr2 * vr2;
      vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      jb_1 = j_1;
      jb_2 = j_2;
      jb_3 = j_3;
      jb_4 = j_4;
      vdqb_1 = vdq_1;
      vdqb_2 = vdq_2;
      vdqb_3 = vdq_3;
      vdqb_4 = vdq_4;
    }
    v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
    v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
    v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
    v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

    v4df vpjb_1 = _mm256_load_pd((double*)(p + jb_1));
    vpjb_1 -= vdf_1 * vdqb_1;
    _mm256_store_pd((double*)(p + jb_1), vpjb_1);

    v4df vpjb_2 = _mm256_load_pd((double*)(p + jb_2));
    vpjb_2 -= vdf_2 * vdqb_2;
    _mm256_store_pd((double*)(p + jb_2), vpjb_2);

    v4df vpjb_3 = _mm256_load_pd((double*)(p + jb_3));
    vpjb_3 -= vdf_3 * vdqb_3;
    _mm256_store_pd((double*)(p + jb_3), vpjb_3);

    v4df vpjb_4 = _mm256_load_pd((double*)(p + jb_4));
    vpjb_4 -= vdf_4 * vdqb_4;
    _mm256_store_pd((double*)(p + jb_4), vpjb_4);

    v4df vpi = _mm256_load_pd((double*)(p + i));
    vpf += vdf_1 * vdqb_1;
    vpf += vdf_2 * vdqb_2;
    vpf += vdf_3 * vdqb_3;
    vpf += vdf_4 * vdqb_4;
    vpi += vpf;
    _mm256_store_pd((double*)(p + i), vpi);
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    double pfx = 0.0;
    double pfy = 0.0;
    double pfz = 0.0;
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[k + kp];
      double dx = q[j][X] - qix;
      double dy = q[j][Y] - qiy;
      double dz = q[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
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
__attribute__((noinline))
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
__attribute__((noinline))
void
force_sorted_swp_intrin_mat_transp(void) {
  const int pn = particle_number;
  const v4df vzero = _mm256_set_pd(0, 0, 0, 0);
  const v4df vcl2 = _mm256_set_pd(CL2, CL2, CL2, CL2);
  const v4df vc24 = _mm256_set_pd(24 * dt, 24 * dt, 24 * dt, 24 * dt);
  const v4df vc48 = _mm256_set_pd(48 * dt, 48 * dt, 48 * dt, 48 * dt);
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(q + i));
    v4df vpf = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    const int kp = pointer[i];
    int ja_1 = sorted_list[kp];
    int ja_2 = sorted_list[kp + 1];
    int ja_3 = sorted_list[kp + 2];
    int ja_4 = sorted_list[kp + 3];
    v4df vqj_1 = _mm256_load_pd((double*)(q + ja_1));
    v4df vdqa_1 = vqj_1 - vqi;
    v4df vqj_2 = _mm256_load_pd((double*)(q + ja_2));
    v4df vdqa_2 = vqj_2 - vqi;
    v4df vqj_3 = _mm256_load_pd((double*)(q + ja_3));
    v4df vdqa_3 = vqj_3 - vqi;
    v4df vqj_4 = _mm256_load_pd((double*)(q + ja_4));
    v4df vdqa_4 = vqj_4 - vqi;

    v4df vdf = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);

    v4df vdqb_1 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_2 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_3 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_4 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);

    int jb_1 = 0, jb_2 = 0, jb_3 = 0, jb_4 = 0;
    const int np = number_of_partners[i];
    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_1 = ja_1;
      const int j_2 = ja_2;
      const int j_3 = ja_3;
      const int j_4 = ja_4;
      v4df vdq_1 = vdqa_1;
      v4df vdq_2 = vdqa_2;
      v4df vdq_3 = vdqa_3;
      v4df vdq_4 = vdqa_4;

      ja_1 = sorted_list[kp + k + 4];
      ja_2 = sorted_list[kp + k + 5];
      ja_3 = sorted_list[kp + k + 6];
      ja_4 = sorted_list[kp + k + 7];

      v4df tmp0 = _mm256_unpacklo_pd(vdq_1, vdq_2);
      v4df tmp1 = _mm256_unpackhi_pd(vdq_1, vdq_2);
      v4df tmp2 = _mm256_unpacklo_pd(vdq_3, vdq_4);
      v4df tmp3 = _mm256_unpackhi_pd(vdq_3, vdq_4);

      v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

      v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

      vqj_1 = _mm256_load_pd((double*)(q + ja_1));
      vdqa_1 = vqj_1 - vqi;
      vpf += vdf_1 * vdqb_1;

      v4df vpjb_1 = _mm256_load_pd((double*)(p + jb_1));
      vpjb_1 -= vdf_1 * vdqb_1;
      _mm256_store_pd((double*)(p + jb_1), vpjb_1);

      vqj_2 = _mm256_load_pd((double*)(q + ja_2));
      vdqa_2 = vqj_2 - vqi;
      vpf += vdf_2 * vdqb_2;

      v4df vpjb_2 = _mm256_load_pd((double*)(p + jb_2));
      vpjb_2 -= vdf_2 * vdqb_2;
      _mm256_store_pd((double*)(p + jb_2), vpjb_2);

      vqj_3 = _mm256_load_pd((double*)(q + ja_3));
      vdqa_3 = vqj_3 - vqi;
      vpf += vdf_3 * vdqb_3;

      v4df vpjb_3 = _mm256_load_pd((double*)(p + jb_3));
      vpjb_3 -= vdf_3 * vdqb_3;
      _mm256_store_pd((double*)(p + jb_3), vpjb_3);

      vqj_4 = _mm256_load_pd((double*)(q + ja_4));
      vdqa_4 = vqj_4 - vqi;
      vpf += vdf_4 * vdqb_4;

      v4df vpjb_4 = _mm256_load_pd((double*)(p + jb_4));
      vpjb_4 -= vdf_4 * vdqb_4;
      _mm256_store_pd((double*)(p + jb_4), vpjb_4);

      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      v4df vr6 = vr2 * vr2 * vr2;
      vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      jb_1 = j_1;
      jb_2 = j_2;
      jb_3 = j_3;
      jb_4 = j_4;
      vdqb_1 = vdq_1;
      vdqb_2 = vdq_2;
      vdqb_3 = vdq_3;
      vdqb_4 = vdq_4;
    }
    v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
    v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
    v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
    v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

    v4df vpjb_1 = _mm256_load_pd((double*)(p + jb_1));
    vpjb_1 -= vdf_1 * vdqb_1;
    _mm256_store_pd((double*)(p + jb_1), vpjb_1);

    v4df vpjb_2 = _mm256_load_pd((double*)(p + jb_2));
    vpjb_2 -= vdf_2 * vdqb_2;
    _mm256_store_pd((double*)(p + jb_2), vpjb_2);

    v4df vpjb_3 = _mm256_load_pd((double*)(p + jb_3));
    vpjb_3 -= vdf_3 * vdqb_3;
    _mm256_store_pd((double*)(p + jb_3), vpjb_3);

    v4df vpjb_4 = _mm256_load_pd((double*)(p + jb_4));
    vpjb_4 -= vdf_4 * vdqb_4;
    _mm256_store_pd((double*)(p + jb_4), vpjb_4);

    v4df vpi = _mm256_load_pd((double*)(p + i));
    vpf += vdf_1 * vdqb_1;
    vpf += vdf_2 * vdqb_2;
    vpf += vdf_3 * vdqb_3;
    vpf += vdf_4 * vdqb_4;
    vpi += vpf;
    _mm256_store_pd((double*)(p + i), vpi);
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    double pfx = 0.0;
    double pfy = 0.0;
    double pfz = 0.0;
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[k + kp];
      double dx = q[j][X] - qix;
      double dy = q[j][Y] - qiy;
      double dz = q[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
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
__attribute__((noinline))
void
measure(void(*pfunc)(), const char *name) {
  double st = myclock();
  const int LOOP = 100;
  // const int LOOP = 1;
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
void
print_result_z(void) {
  for (int i = 0; i < 5; i++) {
    printf("%.10f %.10f %.10f\n", z[i][PX], z[i][PY], z[i][PZ]);
  }
  for (int i = particle_number - 5; i < particle_number; i++) {
    printf("%.10f %.10f %.10f\n", z[i][PX], z[i][PY], z[i][PZ]);
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
#elif P_SWP
  measure(&force_pair_swp, "pair_swp");
  print_result();
#elif P_SWP_INTRIN
  measure(&force_pair_swp_intrin, "pair_swp_intrin");
  print_result();
#elif SORTED
  measure(&force_sorted, "sorted");
  print_result();
#elif S_SWP
  measure(&force_sorted_swp, "sorted_swp");
  print_result();
#elif S_INTRIN
  measure(&force_sorted_intrin, "sorted_intrin");
  print_result();
#elif S_SWP_INTRIN
  measure(&force_sorted_swp_intrin, "sorted_swp_intrin");
  print_result();
#elif S_SWP_INTRIN_M_TRANS
  measure(&force_sorted_swp_intrin_mat_transp, "sorted_swp_intrin_mat_transp");
  print_result();
#elif KNL
  copy_to_z();
  measure(&force_sorted_z_intrin,"sorted_z_intrin");
  copy_from_z();
  print_result();
#elif KNL_SWP
  copy_to_z();
  measure(&force_sorted_z_intrin_swp,"sorted_z_intrin_swp");
  copy_from_z();
  print_result();
#elif KNL_GS
  copy_to_z();
  measure(&force_sorted_z_intrin_gs_swp,"sorted_z_intrin");
  copy_from_z();
  print_result();
#else
  measure(&force_pair, "pair");
  measure(&force_pair_swp, "pair_swp");
  measure(&force_pair_swp_intrin, "pair_swp_intrin");
  measure(&force_sorted, "sorted");
  measure(&force_sorted_swp, "sorted_swp");
  measure(&force_sorted_intrin, "sorted_intrin");
  measure(&force_sorted_swp_intrin, "sorted_swp_intrin");
  measure(&force_sorted_swp_intrin_mat_transp, "sorted_swp_intrin_mat_transp");
  copy_to_2();
  measure(&force_sorted2, "sorted2");
  copy_from_2();
  copy_to_z();
  measure(&force_sorted_z_intrin,"sorted_z_intrin");
  copy_from_z();
  copy_to_z();
  measure(&force_sorted_z_intrin_swp,"sorted_z_intrin_swp");
  copy_from_z();
  copy_to_z();
  measure(&force_sorted_z_intrin_gs_swp,"sorted_z_intrin");
  copy_from_z();
#endif
}
//----------------------------------------------------------------------

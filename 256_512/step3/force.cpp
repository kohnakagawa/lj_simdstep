#include <stdio.h>
#include <immintrin.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
//----------------------------------------------------------------------
const double density = 1.0;
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
sortpair(void){
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
        double x = ix*s;
        double y = iy*s;
        double z = iz*s;
        add_particle(x     ,y   ,z);
        add_particle(x     ,y+hs,z+hs);
        add_particle(x+hs  ,y   ,z+hs);
        add_particle(x+hs  ,y+hs,z);
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
force_pair(void){
  for(int k=0;k<number_of_pairs;k++){
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
force_pair_swp(void){
  int k = 0;
  int i_a = i_particles[k];
  int j_a = j_particles[k];
  double dx_b = q[j_a][X] - q[i_a][X];
  double dy_b = q[j_a][Y] - q[i_a][Y];
  double dz_b = q[j_a][Z] - q[i_a][Z];
  double dx_a,dy_a,dz_a;
  int i_b, j_b;
  double df;
  for(k=1;k<number_of_pairs;k++){
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
    if (r2 > CL2) df=0.0;
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
  if (r2 > CL2) df=0.0;
  p[i_a][X] += df * dx_a;
  p[i_a][Y] += df * dy_a;
  p[i_a][Z] += df * dz_a;
  p[j_a][X] -= df * dx_a;
  p[j_a][Y] -= df * dy_a;
  p[j_a][Z] -= df * dz_a;
}
//----------------------------------------------------------------------
void
force_pair_swp_intrin(void){
  int k = 0;
  int i_a1 = i_particles[k];
  int j_a1 = j_particles[k];
  int i_a2 = i_particles[k+1];
  int j_a2 = j_particles[k+1];
  int i_a3 = i_particles[k+2];
  int j_a3 = j_particles[k+2];
  int i_a4 = i_particles[k+3];
  int j_a4 = j_particles[k+3];
  double dx_b1 = q[j_a1][X] - q[i_a1][X];
  double dy_b1 = q[j_a1][Y] - q[i_a1][Y];
  double dz_b1 = q[j_a1][Z] - q[i_a1][Z];
  double dx_b2 = q[j_a2][X] - q[i_a2][X];
  double dy_b2 = q[j_a2][Y] - q[i_a2][Y];
  double dz_b2 = q[j_a2][Z] - q[i_a2][Z];
  double dx_b3 = q[j_a3][X] - q[i_a3][X];
  double dy_b3 = q[j_a3][Y] - q[i_a3][Y];
  double dz_b3 = q[j_a3][Z] - q[i_a3][Z];
  double dx_b4 = q[j_a4][X] - q[i_a4][X];
  double dy_b4 = q[j_a4][Y] - q[i_a4][Y];
  double dz_b4 = q[j_a4][Z] - q[i_a4][Z];
  double dx_a1,dy_a1,dz_a1;
  double dx_a2,dy_a2,dz_a2;
  double dx_a3,dy_a3,dz_a3;
  double dx_a4,dy_a4,dz_a4;
  int i_b1, j_b1;
  int i_b2, j_b2;
  int i_b3, j_b3;
  int i_b4, j_b4;
  double df_1,df_2,df_3,df_4;
  for(k=4;k<(number_of_pairs)/4*4;k+=4){
    dx_a1 = dx_b1; 
    dy_a1 = dy_b1; 
    dz_a1 = dz_b1; 
    dx_a2 = dx_b2; 
    dy_a2 = dy_b2; 
    dz_a2 = dz_b2; 
    dx_a3 = dx_b3; 
    dy_a3 = dy_b3; 
    dz_a3 = dz_b3; 
    dx_a4 = dx_b4; 
    dy_a4 = dy_b4; 
    dz_a4 = dz_b4; 

    i_b1 = i_particles[k];
    j_b1 = j_particles[k];
    i_b2 = i_particles[k+1];
    j_b2 = j_particles[k+1];
    i_b3 = i_particles[k+2];
    j_b3 = j_particles[k+2];
    i_b4 = i_particles[k+3];
    j_b4 = j_particles[k+3];

    dx_b1 = q[j_b1][X] - q[i_b1][X];
    dy_b1 = q[j_b1][Y] - q[i_b1][Y];
    dz_b1 = q[j_b1][Z] - q[i_b1][Z];
    dx_b2 = q[j_b2][X] - q[i_b2][X];
    dy_b2 = q[j_b2][Y] - q[i_b2][Y];
    dz_b2 = q[j_b2][Z] - q[i_b2][Z];
    dx_b3 = q[j_b3][X] - q[i_b3][X];
    dy_b3 = q[j_b3][Y] - q[i_b3][Y];
    dz_b3 = q[j_b3][Z] - q[i_b3][Z];
    dx_b4 = q[j_b4][X] - q[i_b4][X];
    dy_b4 = q[j_b4][Y] - q[i_b4][Y];
    dz_b4 = q[j_b4][Z] - q[i_b4][Z];

    const double r2_1 = (dx_a1 * dx_a1 + dy_a1 * dy_a1 + dz_a1 * dz_a1);
    const double r2_2 = (dx_a2 * dx_a2 + dy_a2 * dy_a2 + dz_a2 * dz_a2);
    const double r2_3 = (dx_a3 * dx_a3 + dy_a3 * dy_a3 + dz_a3 * dz_a3);
    const double r2_4 = (dx_a4 * dx_a4 + dy_a4 * dy_a4 + dz_a4 * dz_a4);
    const double r6_1 = r2_1 * r2_1 * r2_1;
    const double r6_2 = r2_2 * r2_2 * r2_2;
    const double r6_3 = r2_3 * r2_3 * r2_3;
    const double r6_4 = r2_4 * r2_4 * r2_4;
    df_1 = ((24.0 * r6_1 - 48.0) / (r6_1 * r6_1 * r2_1)) * dt;
    df_2 = ((24.0 * r6_2 - 48.0) / (r6_2 * r6_2 * r2_2)) * dt;
    df_3 = ((24.0 * r6_3 - 48.0) / (r6_3 * r6_3 * r2_3)) * dt;
    df_4 = ((24.0 * r6_4 - 48.0) / (r6_4 * r6_4 * r2_4)) * dt;

    if (r2_1 > CL2) df_1=0.0;
    if (r2_2 > CL2) df_2=0.0;
    if (r2_3 > CL2) df_3=0.0;
    if (r2_4 > CL2) df_4=0.0;
    p[i_a1][X] += df_1 * dx_a1;
    p[i_a1][Y] += df_1 * dy_a1;
    p[i_a1][Z] += df_1 * dz_a1;
    p[j_a1][X] -= df_1 * dx_a1;
    p[j_a1][Y] -= df_1 * dy_a1;
    p[j_a1][Z] -= df_1 * dz_a1;

    p[i_a2][X] += df_2 * dx_a2;
    p[i_a2][Y] += df_2 * dy_a2;
    p[i_a2][Z] += df_2 * dz_a2;
    p[j_a2][X] -= df_2 * dx_a2;
    p[j_a2][Y] -= df_2 * dy_a2;
    p[j_a2][Z] -= df_2 * dz_a2;

    p[i_a3][X] += df_3 * dx_a3;
    p[i_a3][Y] += df_3 * dy_a3;
    p[i_a3][Z] += df_3 * dz_a3;
    p[j_a3][X] -= df_3 * dx_a3;
    p[j_a3][Y] -= df_3 * dy_a3;
    p[j_a3][Z] -= df_3 * dz_a3;

    p[i_a4][X] += df_4 * dx_a4;
    p[i_a4][Y] += df_4 * dy_a4;
    p[i_a4][Z] += df_4 * dz_a4;
    p[j_a4][X] -= df_4 * dx_a4;
    p[j_a4][Y] -= df_4 * dy_a4;
    p[j_a4][Z] -= df_4 * dz_a4;

    i_a1 = i_b1;
    j_a1 = j_b1;
    i_a2 = i_b2;
    j_a2 = j_b2;
    i_a3 = i_b3;
    j_a3 = j_b3;
    i_a4 = i_b4;
    j_a4 = j_b4;
  }
  for(k=(number_of_pairs)/4*4-4;k<number_of_pairs;k++){
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
force_sorted(void){
  const int pn =particle_number;
  for (int i=0; i<pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k=0; k<np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx*dx + dy*dy + dz*dz);
      if (r2 > CL2) continue;
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2))*dt;
      pfx += df*dx;
      pfy += df*dy;
      pfz += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
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
loadpair(void){
  std::ifstream ifs("pair.dat",std::ios::binary);
  ifs.read((char*)&number_of_pairs,sizeof(int));
  ifs.read((char*)number_of_partners,sizeof(int)*N);
  ifs.read((char*)i_particles,sizeof(int)*MAX_PAIRS);
  ifs.read((char*)j_particles,sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
savepair(void){
  makepair();
  std::ofstream ofs("pair.dat",std::ios::binary);
  ofs.write((char*)&number_of_pairs,sizeof(int));
  ofs.write((char*)number_of_partners,sizeof(int)*N);
  ofs.write((char*)i_particles,sizeof(int)*MAX_PAIRS);
  ofs.write((char*)j_particles,sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
print_result(void){
  for (int i = 0; i < 5; i++) {
    printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
  }
  for (int i = particle_number-5; i < particle_number; i++) {
    printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
  }
}
//----------------------------------------------------------------------
int
main(void) {
  init();
  struct stat st;
  int ret = stat("pair.dat",&st);
  if(ret==0){
    std::cerr << "A pair-file is found. I use it." << std::endl;
    loadpair();
  }else{
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
#else
  measure(&force_pair, "pair");
  measure(&force_pair_swp, "pair_swp");
  measure(&force_pair_swp_intrin, "pair_swp_intrin");
  measure(&force_sorted, "sorted");
  measure(&force_sorted_swp, "sorted_swp");
  measure(&force_sorted_intrin, "sorted_intrin");
#endif

}
//----------------------------------------------------------------------

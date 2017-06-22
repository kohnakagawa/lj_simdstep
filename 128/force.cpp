#include <stdio.h>
#include <algorithm>
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
enum {X, Y, Z, W};
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
typedef double v2df __attribute__((vector_size(16)));
//----------------------------------------------------------------------
void
random_shfl() {
  std::mt19937 mt(10);
  const auto pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const auto kp = pointer[i];
    const auto np = number_of_partners[i];
    std::shuffle(&sorted_list[kp], &sorted_list[kp + np], mt);
  }
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
void
force_sorted_intrin_gs(void) {
  const v2df vzero = _mm_set1_pd(0);
  const v2df vcl2  = _mm_set1_pd(CL2);
  const v2df vc24  = _mm_set1_pd(24.0 * dt);
  const v2df vc48  = _mm_set1_pd(48.0 * dt);
  const int pn     = particle_number;
  for (int i = 0; i < pn; i++) {
    const v2df vqxi     = _mm_set1_pd(q[i][X]);
    const v2df vqyi     = _mm_set1_pd(q[i][Y]);
    const v2df vqzi     = _mm_set1_pd(q[i][Z]);
    v2df vpxi           = _mm_setzero_pd();
    v2df vpyi           = _mm_setzero_pd();
    v2df vpzi           = _mm_setzero_pd();
    const int np        = number_of_partners[i];
    const int kp        = pointer[i];
    const int* ptr_list = &sorted_list[kp];
    for (int k = 0; k < (np / 2) * 2; k += 2) {
      const int j_a = *ptr_list++;
      const int j_b = *ptr_list++;
      const v2df vqxj = _mm_set_pd(q[j_b][X], q[j_a][X]);
      const v2df vqyj = _mm_set_pd(q[j_b][Y], q[j_a][Y]);
      const v2df vqzj = _mm_set_pd(q[j_b][Z], q[j_a][Z]);

      _mm256_zeroupper();

      const v2df vdx = vqxj - vqxi;
      const v2df vdy = vqyj - vqyi;
      const v2df vdz = vqzj - vqzi;

      const v2df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      const v2df vr6 = vr2 * vr2 * vr2;
      v2df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      const v2df mask = vcl2 - vr2;
      vdf = _mm_blendv_pd(vdf, vzero, mask);

      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      _mm256_zeroupper();

      v2df vpxj = _mm_set_pd(p[j_b][X], p[j_a][X]);
      v2df vpyj = _mm_set_pd(p[j_b][Y], p[j_a][Y]);
      v2df vpzj = _mm_set_pd(p[j_b][Z], p[j_a][Z]);

      _mm256_zeroupper();

      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;

      const double pxj_a = _mm_cvtsd_f64(vpxj);
      const double pyj_a = _mm_cvtsd_f64(vpyj);
      const double pzj_a = _mm_cvtsd_f64(vpzj);

      const double pxj_b = _mm_cvtsd_f64(_mm_castsi128_pd(_mm_bsrli_si128(_mm_castpd_si128(vpxj), 8)));
      const double pyj_b = _mm_cvtsd_f64(_mm_castsi128_pd(_mm_bsrli_si128(_mm_castpd_si128(vpyj), 8)));
      const double pzj_b = _mm_cvtsd_f64(_mm_castsi128_pd(_mm_bsrli_si128(_mm_castpd_si128(vpzj), 8)));

      p[j_a][X] = pxj_a;
      p[j_a][Y] = pyj_a;
      p[j_a][Z] = pzj_a;

      p[j_b][X] = pxj_b;
      p[j_b][Y] = pyj_b;
      p[j_b][Z] = pzj_b;
    }
    auto pix = p[i][X];
    auto piy = p[i][Y];
    auto piz = p[i][Z];

    pix += _mm_cvtsd_f64(vpxi);
    piy += _mm_cvtsd_f64(vpyi);
    piz += _mm_cvtsd_f64(vpzi);

    vpxi = _mm_castsi128_pd(_mm_bsrli_si128(_mm_castpd_si128(vpxi), 8));
    vpyi = _mm_castsi128_pd(_mm_bsrli_si128(_mm_castpd_si128(vpyi), 8));
    vpzi = _mm_castsi128_pd(_mm_bsrli_si128(_mm_castpd_si128(vpzi), 8));

    pix += _mm_cvtsd_f64(vpxi);
    piy += _mm_cvtsd_f64(vpyi);
    piz += _mm_cvtsd_f64(vpzi);

    p[i][X] = pix;
    p[i][Y] = piy;
    p[i][Z] = piz;
    if (np & 0x1) {
      const int j = *ptr_list;
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
force_sorted_intrin_nogs(void) {
  const v2df vzero = _mm_set1_pd(0);
  const v2df vcl2  = _mm_set1_pd(CL2);
  const v2df vc24  = _mm_set1_pd(24.0 * dt);
  const v2df vc48  = _mm_set1_pd(48.0 * dt);
  const int pn     = particle_number;
  for (int i = 0; i < pn; i++) {
    const v2df vqi_yx   = _mm_load_pd(&q[i][X]);
    const v2df vqi_wz   = _mm_load_pd(&q[i][Z]);
    v2df vpi_yx         = _mm_setzero_pd();
    v2df vpi_wz         = _mm_setzero_pd();
    const int np        = number_of_partners[i];
    const int kp        = pointer[i];
    const int* ptr_list = &sorted_list[kp];
    for (int k = 0; k < (np / 2) * 2; k += 2) {
      const int j_a = *ptr_list++;
      const int j_b = *ptr_list++;

      v2df vqja_yx = _mm_load_pd(&q[j_a][X]);
      v2df vqja_wz = _mm_load_pd(&q[j_a][Z]);
      v2df vqjb_yx = _mm_load_pd(&q[j_b][X]);
      v2df vqjb_wz = _mm_load_pd(&q[j_b][Z]);

      v2df vdqa_yx = vqja_yx - vqi_yx;
      v2df vdqa_wz = vqja_wz - vqi_wz;
      v2df vdqb_yx = vqjb_yx - vqi_yx;
      v2df vdqb_wz = vqjb_wz - vqi_wz;

      v2df vdx = _mm_unpacklo_pd(vdqa_yx, vdqb_yx);
      v2df vdy = _mm_unpackhi_pd(vdqa_yx, vdqb_yx);
      v2df vdz = _mm_unpacklo_pd(vdqa_wz, vdqb_wz);

      v2df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      v2df vr6 = vr2 * vr2 * vr2;
      v2df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v2df mask = vcl2 - vr2;
      vdf = _mm_blendv_pd(vdf, vzero, mask);

      v2df vdf_a = _mm_permute_pd(vdf, 0x0);
      v2df vdf_b = _mm_permute_pd(vdf, 0x3);

      vpi_yx += vdf_a * vdqa_yx;
      vpi_yx += vdf_b * vdqb_yx;
      vpi_wz += vdf_a * vdqa_wz;
      vpi_wz += vdf_b * vdqb_wz;

      v2df vpja_yx = _mm_load_pd(&p[j_a][X]);
      v2df vpja_wz = _mm_load_pd(&p[j_a][Z]);
      v2df vpjb_yx = _mm_load_pd(&p[j_b][X]);
      v2df vpjb_wz = _mm_load_pd(&p[j_b][Z]);

      vpja_yx -= vdf_a * vdqa_yx;
      vpjb_yx -= vdf_b * vdqb_yx;
      vpja_wz -= vdf_a * vdqa_wz;
      vpjb_wz -= vdf_b * vdqb_wz;

      _mm_store_pd(&p[j_a][X], vpja_yx);
      _mm_store_pd(&p[j_a][Z], vpja_wz);
      _mm_store_pd(&p[j_b][X], vpjb_yx);
      _mm_store_pd(&p[j_b][Z], vpjb_wz);
    }
    vpi_yx += (v2df)_mm_load_pd(&p[i][X]);
    vpi_wz += (v2df)_mm_load_pd(&p[i][Z]);
    _mm_store_pd(&p[i][X], vpi_yx);
    _mm_store_pd(&p[i][Z], vpi_wz);

    if (np & 0x1) {
      const int j = *ptr_list;
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
  //const int LOOP = 1;
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
#ifdef SORTED
  measure(&force_sorted, "sorted");
  print_result();
#elif S_INTRIN_GS
  measure(&force_sorted_intrin_gs, "sorted_intrin_gs");
  print_result();
#elif S_INTRIN_NOGS
  measure(&force_sorted_intrin_nogs, "sorted_intrin_nogs");
  print_result();
#else
  measure(&force_sorted, "sorted");
  measure(&force_sorted_intrin_gs, "sorted_intrin_gs");
  measure(&force_sorted_intrin_nogs, "sorted_intrin_nogs");
#endif
}
//----------------------------------------------------------------------

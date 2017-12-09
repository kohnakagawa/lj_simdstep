#include <iostream>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <random>
#include <cmath>
#include <chrono>
#include <x86intrin.h>
#include <sys/stat.h>
//----------------------------------------------------------------------
struct double4 {double x, y, z, w;};
struct double8 {double x, y, z, w, px, py, pz, pw;};
enum {X = 0, Y, Z, W, PX, PY, PZ, PW, DIM};
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_plain(const double4* q,
            double4* p,
            const int* number_of_partners,
            const int* pointer,
            const int* sorted_list,
            const double CL2,
            const double dt,
            const int pn) {
  for (int i = 0; i < pn; i++) {
    const auto qx_key = q[i].x;
    const auto qy_key = q[i].y;
    const auto qz_key = q[i].z;
    const auto np = number_of_partners[i];
    double pfx = 0, pfy = 0, pfz = 0;
    const auto kp = pointer[i];
#pragma novector
    for (int k = 0; k < np; k++) {
      const auto j = sorted_list[kp + k];
      const auto dx = q[j].x - qx_key;
      const auto dy = q[j].y - qy_key;
      const auto dz = q[j].z - qz_key;
      const auto r2 = (dx*dx + dy*dy + dz*dz);
      if (r2 > CL2) continue;
      const auto r6 = r2*r2*r2;
      const auto df = ((24.0 * r6 - 48.0)/(r6 * r6 * r2)) * dt;
      pfx += df*dx;
      pfy += df*dy;
      pfz += df*dz;
      p[j].x -= df*dx;
      p[j].y -= df*dy;
      p[j].z -= df*dz;
    } // end of k loop
    p[i].x += pfx;
    p[i].y += pfy;
    p[i].z += pfz;
  } // end of i loop
}
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_sorted(const double4* __restrict q,
             double4* __restrict p,
             const int* __restrict number_of_partners,
             const int* __restrict pointer,
             const int* __restrict sorted_list,
             const double CL2,
             const double dt,
             const int pn) {
  for (int i = 0; i < pn; i++) {
    const auto qx_key = q[i].x;
    const auto qy_key = q[i].y;
    const auto qz_key = q[i].z;
    const auto np = number_of_partners[i];
    double pfx = 0, pfy = 0, pfz = 0;
    const auto kp = pointer[i];
#pragma simd
#pragma vector aligned
    for (int k = 0; k < np; k++) {
      const auto j = sorted_list[kp + k];
      const auto dx = q[j].x - qx_key;
      const auto dy = q[j].y - qy_key;
      const auto dz = q[j].z - qz_key;
      const auto r2 = (dx*dx + dy*dy + dz*dz);
      if (r2 > CL2) continue;
      const auto r6 = r2*r2*r2;
      const auto df = ((24.0 * r6 - 48.0)/(r6 * r6 * r2)) * dt;
      pfx += df*dx;
      pfy += df*dy;
      pfz += df*dz;
      p[j].x -= df*dx;
      p[j].y -= df*dy;
      p[j].z -= df*dz;
    } // end of k loop
    p[i].x += pfx;
    p[i].y += pfy;
    p[i].z += pfz;
  } // end of i loop
}
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_sorted_z(double8* __restrict z,
               const int* __restrict number_of_partners,
               const int* __restrict pointer,
               const int* __restrict sorted_list,
               const double CL2,
               const double dt,
               const int pn) {
  for (int i = 0; i < pn; i++) {
    const auto qx_key = z[i].x;
    const auto qy_key = z[i].y;
    const auto qz_key = z[i].z;
    const auto np = number_of_partners[i];
    double pfx = 0, pfy = 0, pfz = 0;
    const auto kp = pointer[i];
#pragma simd
#pragma vector aligned
    for (int k = 0; k < np; k++) {
      const auto j = sorted_list[kp + k];
      const auto dx = z[j].x - qx_key;
      const auto dy = z[j].y - qy_key;
      const auto dz = z[j].z - qz_key;
      const auto r2 = (dx*dx + dy*dy + dz*dz);
      if (r2 > CL2) continue;
      const auto r6 = r2*r2*r2;
      const auto df = ((24.0 * r6 - 48.0)/(r6 * r6 * r2)) * dt;
      pfx += df*dx;
      pfy += df*dy;
      pfz += df*dz;
      z[j].px -= df*dx;
      z[j].py -= df*dy;
      z[j].pz -= df*dz;
    } // end of k loop
    z[i].px += pfx;
    z[i].py += pfy;
    z[i].pz += pfz;
  } // end of i loop
}
//----------------------------------------------------------------------
__attribute__((noinline))
void
force_sorted_z_2d(double z[][8],
                  const int* __restrict number_of_partners,
                  const int* __restrict pointer,
                  const int* __restrict sorted_list,
                  const double CL2,
                  const double dt,
                  const int pn) {
  for (int i = 0; i < pn; i++) {
    const auto qx_key = z[i][X];
    const auto qy_key = z[i][Y];
    const auto qz_key = z[i][Z];
    const auto np = number_of_partners[i];
    double pfx = 0, pfy = 0, pfz = 0;
    const auto kp = pointer[i];
#pragma simd
#pragma vector aligned
    for (int k = 0; k < np; k++) {
      const auto j = sorted_list[kp + k];
      const auto dx = z[j][X] - qx_key;
      const auto dy = z[j][Y] - qy_key;
      const auto dz = z[j][Z] - qz_key;
      const auto r2 = (dx*dx + dy*dy + dz*dz);
      if (r2 > CL2) continue;
      const auto r6 = r2*r2*r2;
      const auto df = ((24.0 * r6 - 48.0)/(r6 * r6 * r2)) * dt;
      pfx += df*dx;
      pfy += df*dy;
      pfz += df*dz;
      z[j][PX] -= df*dx;
      z[j][PY] -= df*dy;
      z[j][PZ] -= df*dz;
    } // end of k loop
    z[i][PX] += pfx;
    z[i][PY] += pfy;
    z[i][PZ] += pfz;
  } // end of i loop
}
//----------------------------------------------------------------------
const double density = 1.0;
//const double density = 0.5;
const int N = 400000;
const int MAX_PAIRS = 30 * N;
const double L = 50.0;
const double dt = 0.001;

const char* pairlist_cache_file_name = "pair.dat";

double4* q = nullptr;
double4* p = nullptr;
double8* z = nullptr;
double z_2d[N][DIM];

int particle_number = 0;
int number_of_pairs = 0;
int* number_of_partners = nullptr;
int i_particles[MAX_PAIRS];
int j_particles[MAX_PAIRS];
int32_t* pointer = nullptr;
int32_t pointer2[N];
int* sorted_list = nullptr;

const double CUTOFF_LENGTH = 3.0;
const double SEARCH_LENGTH = 3.3;
const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
//----------------------------------------------------------------------
void
add_particle(double x, double y, double z) {
  static std::mt19937 mt(2);
  std::uniform_real_distribution<double> ud(0.0, 0.1);
  q[particle_number].x = x + ud(mt);
  q[particle_number].y = y + ud(mt);
  q[particle_number].z = z + ud(mt);
  particle_number++;
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
      const double dx = q[i].x - q[j].x;
      const double dy = q[i].y - q[j].y;
      const double dz = q[i].z - q[j].z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < SL2) {
        register_pair(i, j);
      }
    }
  }
}
//----------------------------------------------------------------------
void
allocate(void) {
  posix_memalign((void**)(&q), 64, sizeof(double4) * N);
  posix_memalign((void**)(&p), 64, sizeof(double4) * N);
  posix_memalign((void**)(&z), 64, sizeof(double8) * N);

  posix_memalign((void**)(&number_of_partners), 64, sizeof(int) * N);
  posix_memalign((void**)(&pointer), 64, sizeof(int32_t) * N);
  posix_memalign((void**)(&sorted_list), 64, sizeof(int32_t) * MAX_PAIRS);

  std::fill(number_of_partners,
            number_of_partners + N,
            0);
  std::fill(pointer,
            pointer + N,
            0);
  std::fill(sorted_list,
            sorted_list + MAX_PAIRS,
            0);
}
//----------------------------------------------------------------------
void
deallocate(void) {
  free(q);
  free(p);
  free(z);

  free(number_of_partners);
  free(pointer);
  free(sorted_list);
}
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
copy(void) {
  for (int i = 0; i < particle_number; i++) {
    z_2d[i][X] = z[i].x = q[i].x;
    z_2d[i][Y] = z[i].y = q[i].y;
    z_2d[i][Z] = z[i].z = q[i].z;
    z_2d[i][W] = z[i].w = 0.0;

    z_2d[i][PX] = z[i].px = p[i].x;
    z_2d[i][PY] = z[i].py = p[i].y;
    z_2d[i][PZ] = z[i].pz = p[i].z;
    z_2d[i][PW] = z[i].pw = 0.0;
  }
}
//----------------------------------------------------------------------
void
init(void) {
  const double s = 1.0 / std::pow(density * 0.25, 1.0 / 3.0);
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
    p[i].x = 0.0;
    p[i].y = 0.0;
    p[i].z = 0.0;
  }
}
//----------------------------------------------------------------------
#define MEASURE(func, name)						\
  do {									\
    const auto beg = std::chrono::system_clock::now();			\
    const int LOOP = 100;						\
    for (int i = 0; i < LOOP; i++) {					\
      func;								\
    }									\
    const auto end = std::chrono::system_clock::now();			\
    const long dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg).count(); \
    fprintf(stderr, "N=%d, %s %ld [ms]\n", particle_number, name, dur);	\
  } while (0)
//----------------------------------------------------------------------
void
loadpair(void){
  std::ifstream ifs(pairlist_cache_file_name,std::ios::binary);
  ifs.read((char*)&number_of_pairs,sizeof(int));
  ifs.read((char*)number_of_partners,sizeof(int)*N);
  ifs.read((char*)i_particles,sizeof(int)*MAX_PAIRS);
  ifs.read((char*)j_particles,sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
savepair(void){
  makepair();
  random_shfl();
  std::ofstream ofs(pairlist_cache_file_name,std::ios::binary);
  ofs.write((char*)&number_of_pairs,sizeof(int));
  ofs.write((char*)number_of_partners,sizeof(int)*N);
  ofs.write((char*)i_particles,sizeof(int)*MAX_PAIRS);
  ofs.write((char*)j_particles,sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
print_result_4(double4* mom){
  for (int i = 0; i < 5; i++) {
    printf("%.10f %.10f %.10f\n", mom[i].x, mom[i].y, mom[i].z);
  }
  for (int i = particle_number-5; i < particle_number; i++) {
    printf("%.10f %.10f %.10f\n", mom[i].x, mom[i].y, mom[i].z);
  }
  printf("\n");
}
//----------------------------------------------------------------------
void
print_result_8(double8* mom) {
  for (int i = 0; i < 5; i++) {
    printf("%.10f %.10f %.10f\n", mom[i].px, mom[i].py, mom[i].pz);
  }
  for (int i = particle_number-5; i < particle_number; i++) {
    printf("%.10f %.10f %.10f\n", mom[i].px, mom[i].py, mom[i].pz);
  }
  printf("\n");
}

//----------------------------------------------------------------------
int
main(void) {
  allocate();
  init();
  struct stat st;
  int ret = stat(pairlist_cache_file_name, &st);
  if (ret == 0) {
    std::cerr << "A pair-file is found. I use it." << std::endl;
    loadpair();
  } else {
    std::cerr << "Make pairlist." << std::endl;
    savepair();
  }
  std::cerr << "Number of pairs: " << number_of_pairs << std::endl;

  sortpair();
  copy();

  MEASURE(force_plain(q, p, number_of_partners, pointer, sorted_list, CL2, dt, particle_number), "plain");
  print_result_4(p);

  MEASURE(force_sorted(q, p, number_of_partners, pointer, sorted_list, CL2, dt, particle_number), "sorted");
  print_result_4(p);

  MEASURE(force_sorted_z(z, number_of_partners, pointer, sorted_list, CL2, dt, particle_number), "sorted_z");
  print_result_8(z);

  MEASURE(force_sorted_z_2d(z_2d, number_of_partners, pointer, sorted_list, CL2, dt, particle_number), "sorted_z_2d");
  print_result_8((double8*)&z_2d[0][0]);

  deallocate();
}
//----------------------------------------------------------------------

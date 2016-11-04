#include <stdio.h>
#include <immintrin.h>
//----------------------------------------------------------------------
typedef double v4df __attribute__((vector_size(32)));
//----------------------------------------------------------------------
#define p4(x) printf("%f %f\n",x##1,x##2)
//----------------------------------------------------------------------
void
print256(v4df r) {
  double *a = (double*)(&r);
  printf("%.10f %.10f %.10f %.10f\n", a[0], a[1], a[2], a[3]);
}
//----------------------------------------------------------------------
int
main(void){
  double x[] = {0,1,2,3};
  v4df r1 = _mm256_load_pd(x);
  v4df r2 = _mm256_set_pd(x[0],x[1],x[2],x[3]);
  print256(r1);
  print256(r2);
  v4df r3 = _mm256_permute4x64_pd(r1, 0);
  print256(r3);
  v4df r4 = _mm256_permute4x64_pd(r2, 0);
  print256(r4);
  double y1 = 1.0;
  double y2 = 2.0;
  p4(y);
}

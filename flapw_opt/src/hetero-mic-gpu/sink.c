#include <hStreams_sink.h>
#include <stdio.h>
#include <time.h>
#include <mkl.h>

typedef struct cx_double {
  double re, im;
} cx_double;

void zher2k_(const char *, const char *, const int *, const int *, const cx_double *, const cx_double *, const int *, const cx_double *, const int *, const cx_double *, cx_double *, const int *);
void zherk_(const char *, const char *, const int *, const int *, const cx_double *, const cx_double *, const int *, const cx_double *, cx_double *, const int *);
void zgemm_(const char *, const char *, const int *, const int *, const int *, const cx_double *, const cx_double *, const int *, const cx_double *, const int *, const cx_double *, cx_double *, const int *);

static const cx_double c0 = {0, 0};
static const cx_double c1 = {1, 0};

#define HS_GEMM 0
#define HS_HER2K 1
#define HS_HERK 2
#define HS_HERKX 3

HSTREAMS_EXPORT
void hs_thunk_sink(
  uint64_t arg0, uint64_t arg1, uint64_t arg2, uint64_t arg3, uint64_t arg4, uint64_t arg5, uint64_t arg6, uint64_t arg7
) {
  int kind       =               arg0;
  int dev        =               arg1;
  int theM       =               arg2;
  int theN       =               arg3;
  int theK       =               arg4;
  cx_double * bA = (cx_double *) arg5;
  cx_double * bB = (cx_double *) arg6;
  cx_double * bC = (cx_double *) arg7;
  struct timespec ts, te;
  clock_gettime(CLOCK_REALTIME, &ts);
  double start = ts.tv_sec + 1.0e-9 * ts.tv_nsec;
  if (kind == HS_GEMM) {
    zgemm_("C", "N", &theM, &theN, &theK, &c1, bA, &theK, bB, &theK, &c0, bC, &theM);
  } else if (kind == HS_HER2K) {
    zher2k_("L", "C", &theN, &theK, &c1, bA, &theK, bB, &theK, &c0, bC, &theN);
  } else if (kind == HS_HERK) {
    zherk_("L", "C", &theN, &theK, &c1, bA, &theK, &c0, bC, &theN);
  } else if (kind == HS_HERKX) {
    zgemm_("C", "N", &theM, &theN, &theK, &c1, bA, &theK, bB, &theK, &c0, bC, &theM);
    //zgemmt_("L", "C", "N", &theN, &theK, &c1, bA, &theK, bB, &theK, &c0, bC, &theN);
    int i;
    for (i = 0; i < theN; i++) { // prob not worth parallelizing.
      bC[i * theM + i].im = 0;
    }
  }
  clock_gettime(CLOCK_REALTIME, &te);
  double end = te.tv_sec + 1.0e-9 * te.tv_nsec;
  double elapsed = te.tv_sec - ts.tv_sec + 1.0e-9 * (te.tv_nsec - ts.tv_nsec);
  fprintf(stderr, "THREADS %d\n", mkl_get_max_threads());
  fprintf(stderr, "BLOCK FLOPS %e\n", (kind == HS_HERK ? 4.0 : 8.0) * theN * theM * theK / elapsed);
  fprintf(stderr, "TRACE time dev%d %.20e %.20e\n", dev, start, end);
}

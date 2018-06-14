#if defined(HS_CUDA) || defined(HS_HSTREAMS)
#include <assert.h>
#include <stdint.h>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <queue>
#include <memory>
#include <omp.h>
#ifdef HS_HSTREAMS
# include <hStreams_source.h>
# include <hStreams_common.h>
# include <hStreams_app_api.h>
#elif defined(HS_CUDA)
# include <cuda_runtime.h>
# include <cublas_v2.h>
#else
# error "Need to compile with either hetero streams or cuda"
#endif

#include "blas.h"
#include "source.h"

#define LOG_HSTR_RESULT(func)                                                 \
    {                                                                         \
        HSTR_RESULT hret = HSTR_RESULT_SUCCESS;                               \
        hret = func;                                                          \
        if (hret != HSTR_RESULT_SUCCESS) {                                    \
            printf("%s returned %s. (%s:%d)\n", #func, hStreams_ResultGetName(hret), __FILE__, __LINE__); \
            abort();                                                           \
        }                                                                     \
    }

#define LOG_CUDA_ERROR(func)                                                 \
    {                                                                         \
        cudaError_t hret = cudaSuccess;                               \
        hret = func;                                                          \
        if (hret != cudaSuccess) {                                    \
            printf("%s returned %s. (%s:%d)\n", #func, cudaGetErrorString(hret), __FILE__, __LINE__); \
            abort();                                                           \
        }                                                                     \
    }

#ifdef HS_CUDA
const char * cublasGetErrorString(cublasStatus_t status) {
  // Stolen: http://stackoverflow.com/questions/13041399/equivalent-of-cudageterrorstring-for-cublas
  switch(status) {
    case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE"; 
    case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH"; 
    case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED"; 
    case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR"; 
  }
  return "unknown error";
}
#endif

#define LOG_CUBLAS_STATUS(func)                                                 \
    {                                                                         \
        cublasStatus_t hret = CUBLAS_STATUS_SUCCESS;                               \
        hret = func;                                                          \
        if (hret != CUBLAS_STATUS_SUCCESS) {                                    \
            printf("%s returned %s. (%s:%d)\n", #func, cublasGetErrorString(hret), __FILE__, __LINE__); \
            abort();                                                           \
        }                                                                     \
    }

#define HS_GEMM 0
#define HS_HER2K 1
#define HS_HERK 2
#define HS_HERKX 3

struct hs_clock {
  double last;
  double first;
  int has_first;

  double the_time() {
    struct timespec ret;
    clock_gettime(CLOCK_REALTIME, &ret);
    return ret.tv_sec + 1.0e-9 * ret.tv_nsec;
  }

  hs_clock() : last(the_time()), has_first(0) {}

  double mark() {
    double ret = the_time();
    if (! has_first) {
      first = ret;
       has_first = 1;
    }
    assert(ret >= last);
    last = ret;
    return ret;
  }
} hs_g_clock;

struct hs_timing_data {
  double sum;
  std::vector<double> individual;
  std::vector<double> start;
  std::vector<double> end;
  int len;

  hs_timing_data() : sum(0), len(0) {}

  void output(const char * name) {
    if (sum >= 1)
      printf("%s sum %5.1f\n", name, sum);
    if (getenv("HSTRACE") && strcmp("1", getenv("HSTRACE")) == 0)
      for (int i = 0; i < len; i++) {
        printf("TRACE %s %.20e %.20e\n", name, start[i], end[i]);
      }
  }
};

struct hs_timing_guard {
  struct timespec ts, te;
  double start;
  hs_timing_data &data;

  hs_timing_guard(hs_timing_data &d) : data(d) {
    start = hs_g_clock.mark();
  }

  ~hs_timing_guard() {
    double end = hs_g_clock.mark();
    double diff = end - start;
    data.len += 1;
    data.sum += diff;
    data.individual.push_back(diff);
    data.start.push_back(start);
    data.end.push_back(end);
  }
};

hs_timing_data hs_wait_xfer;
hs_timing_data hs_wait_compute;
hs_timing_data hs_wait_issue_xfer;
hs_timing_data hs_wait_issue_compute;
hs_timing_data hs_wait_issue_event;
hs_timing_data hs_wait_issue_output;
hs_timing_data hs_time_pack;
hs_timing_data hs_time_unpack;
hs_timing_data hs_time_zero;
hs_timing_data hs_check_compute;
hs_timing_data hs_host_compute;
hs_timing_data hs_time_total;

struct hs_accelerator {
  const int device;
  const int ndevs;
# ifdef HS_HSTREAMS
  std::map<cx_double*, HSTR_EVENT> signals;
# else
  std::map<cx_double*, cx_double*> host_to_device;
  std::map<cx_double*, cudaEvent_t> signals;
# endif

  static const int EXEC = 0;
  static const int H2D = 1;
  static const int D2H = 2;

# ifdef HS_CUDA
  cudaStream_t exec, h2d[2], d2h;
  int h2d_idx;
  cublasHandle_t cublas_handle;
# endif

  hs_accelerator(int dev, int ndev) : device(dev), ndevs(ndev) {
#   ifdef HS_CUDA
    LOG_CUDA_ERROR(cudaSetDevice(device));
    LOG_CUDA_ERROR(cudaStreamCreate(&exec));
    LOG_CUDA_ERROR(cudaStreamCreate(&h2d[0]));
    LOG_CUDA_ERROR(cudaStreamCreate(&h2d[1]));
    h2d_idx = 0;
    LOG_CUDA_ERROR(cudaStreamCreate(&d2h));
    LOG_CUBLAS_STATUS(cublasCreate(&cublas_handle));
#   endif
  }

  ~hs_accelerator() {
#   ifdef HS_CUDA
    LOG_CUDA_ERROR(cudaSetDevice(device));
    LOG_CUDA_ERROR(cudaStreamDestroy(exec));
    LOG_CUDA_ERROR(cudaStreamDestroy(h2d[0]));
    LOG_CUDA_ERROR(cudaStreamDestroy(h2d[1]));
    LOG_CUDA_ERROR(cudaStreamDestroy(d2h));
    LOG_CUBLAS_STATUS(cublasDestroy(cublas_handle));
#   endif
  }

  int get_stream_id(int kind) {
    return device + ndevs * kind;
  }

  void wait(cx_double * buf) {
    hs_timing_guard tg(hs_wait_compute);
#   ifdef HS_HSTREAMS
    LOG_HSTR_RESULT(hStreams_app_event_wait(1, &signals[buf]));
#   else
    LOG_CUDA_ERROR(cudaSetDevice(device));
    LOG_CUDA_ERROR(cudaEventSynchronize(signals[buf]));
#   endif
  }

  bool check(cx_double * buf) {
    hs_timing_guard tg(hs_check_compute);
#   ifdef HS_HSTREAMS
    HSTR_RESULT status = hStreams_EventWait(1, &signals[buf], true, 0, nullptr, nullptr);
    if (status == HSTR_RESULT_SUCCESS) {
      return true;
    } else if (status == HSTR_RESULT_TIME_OUT_REACHED) {
      return false;
    } else {
      LOG_HSTR_RESULT(status);
      return false; // never reached
    }
#   else
    LOG_CUDA_ERROR(cudaSetDevice(device));
    cudaError_t status = cudaEventQuery(signals[buf]);
    if (status == cudaSuccess) {
      return true;
    } else if (status == cudaErrorNotReady) {
      return false;
    } else {
      LOG_CUDA_ERROR(status);
      return false; // never reached
    }
#   endif
  }

  void transfer(cx_double * my_buf, int num) {
    hs_timing_guard tg(hs_wait_issue_xfer);
#   ifdef HS_HSTREAMS
    LOG_HSTR_RESULT(hStreams_app_xfer_memory(get_stream_id(H2D), my_buf, my_buf, num * sizeof(cx_double), HSTR_SRC_TO_SINK, &signals[my_buf]));
#   else
    LOG_CUDA_ERROR(cudaSetDevice(device));
    LOG_CUDA_ERROR(cudaMemcpyAsync(host_to_device[my_buf], my_buf, num * sizeof(cx_double), cudaMemcpyHostToDevice, h2d[h2d_idx]));
    LOG_CUDA_ERROR(cudaEventRecord(signals[my_buf], h2d[h2d_idx]));
    h2d_idx = 1 - h2d_idx;
#   endif
  }

  const cx_double * transfer_strided(cx_double * buffer, const cx_double * data, int ld, int num_strides, int num_in_stride, bool try_avoid_buffer) {
#   ifdef HS_HSTREAMS
    try_avoid_buffer = false;
#   endif
    if (! try_avoid_buffer) {
      {
        hs_timing_guard tg(hs_time_pack);
        #pragma omp parallel for
        for (int i = 0; i < num_strides; i++) {
          //for (int j = 0; j < num_in_stride; j++) {
          //  buffer[i * num_in_stride + j] = data[i * ld + j];
          //}
          memcpy(&buffer[i * num_in_stride], &data[i * ld], num_in_stride * sizeof(cx_double));
        }
      }
      transfer(buffer, num_strides*num_in_stride);
    } else {
#   ifdef HS_CUDA
    LOG_CUDA_ERROR(cudaSetDevice(device));
    LOG_CUBLAS_STATUS(cublasSetMatrixAsync(num_in_stride, num_strides, sizeof(cx_double), data, ld, host_to_device[buffer], num_in_stride, h2d[h2d_idx]));
    LOG_CUDA_ERROR(cudaEventRecord(signals[buffer], h2d[h2d_idx]));
    h2d_idx = 1 - h2d_idx;
#   endif
    }
    return buffer;
  }
  
  void thunk(int kind, int theM, int theN, int theK, cx_double * bA, cx_double * bB, cx_double * bC) {
#   ifdef HS_CUDA
    LOG_CUDA_ERROR(cudaSetDevice(device));
#   endif
    {
      hs_timing_guard tg(hs_wait_xfer);
#     ifdef HS_HSTREAMS
      LOG_HSTR_RESULT(hStreams_app_event_wait(1, &signals[bA]));
      LOG_HSTR_RESULT(hStreams_app_event_wait(1, &signals[bB]));
#     else
      LOG_CUDA_ERROR(cudaEventSynchronize(signals[bA]));
      LOG_CUDA_ERROR(cudaEventSynchronize(signals[bB]));
#     endif
    }
    {
      hs_timing_guard tg(hs_wait_issue_compute);
#     ifdef HS_HSTREAMS
      uint64_t args[8];
      args[0] = kind;
      args[1] = device;
      args[2] = theM;
      args[3] = theN;
      args[4] = theK;
      args[5] = (uint64_t) bA;
      args[6] = (uint64_t) bB;
      args[7] = (uint64_t) bC;
      HSTR_EVENT evt;
      LOG_HSTR_RESULT(hStreams_app_invoke(get_stream_id(EXEC), "hs_thunk_sink", 5, 3, args, &evt, NULL, 0));
#     else
      LOG_CUBLAS_STATUS(cublasSetStream(cublas_handle, exec));
      cuDoubleComplex c0 = make_cuDoubleComplex(0, 0);
      cuDoubleComplex c1 = make_cuDoubleComplex(1, 0);
      double d0 = 0;
      double d1 = 1;
      cuDoubleComplex * dA = (cuDoubleComplex *) host_to_device[bA];
      cuDoubleComplex * dB = (cuDoubleComplex *) host_to_device[bB];
      cuDoubleComplex * dC = (cuDoubleComplex *) host_to_device[bC];
      if (kind == HS_HER2K) {
        LOG_CUBLAS_STATUS(cublasZher2k(cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_C, theN, theK, &c1, dA, theK, dB, theK, &d0, dC, theN));
      } else if (kind == HS_HERK) {
        LOG_CUBLAS_STATUS(cublasZherk(cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_C, theN, theK, &d1, dA, theK, &d0, dC, theN));
      } else if (kind == HS_GEMM) {
        LOG_CUBLAS_STATUS(cublasZgemm(cublas_handle, CUBLAS_OP_C, CUBLAS_OP_N, theM, theN, theK, &c1, dA, theK, dB, theK, &c0, dC, theM));
      } else if (kind == HS_HERKX) {
        LOG_CUBLAS_STATUS(cublasZherkx(cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_C, theN, theK, &c1, dA, theK, dB, theK, &d0, dC, theN));
      }
#     endif
    }
    {
      hs_timing_guard tg(hs_wait_issue_output);
#     ifdef HS_HSTREAMS
      LOG_HSTR_RESULT(hStreams_app_xfer_memory(get_stream_id(EXEC), bC, bC, theM*theN*sizeof(cx_double), HSTR_SINK_TO_SRC, &signals[bC]));
      //LOG_HSTR_RESULT(hStreams_EventStreamWait(get_stream_id(device, D2H), 1, &evt, 0, NULL, NULL));
#     else
      LOG_CUDA_ERROR(cudaEventRecord(signals[bC], exec));
      LOG_CUDA_ERROR(cudaStreamWaitEvent(d2h, signals[bC], 0));
      LOG_CUDA_ERROR(cudaMemcpyAsync(bC, host_to_device[bC], theM*theN*sizeof(cx_double), cudaMemcpyDeviceToHost, d2h));
      LOG_CUDA_ERROR(cudaEventRecord(signals[bC], d2h));
      //LOG_CUDA_ERROR(cudaMemcpyAsync(bC, host_to_device[bC], theM*theN*sizeof(cx_double), cudaMemcpyDeviceToHost, exec));
      //LOG_CUDA_ERROR(cudaEventRecord(signals[bC], exec));
#     endif
    }
  }

  cx_double * alloc(int num) {
#   ifdef HS_HSTREAMS
    cx_double * ptr = (cx_double *) malloc(num * sizeof(cx_double));
    size_t partitions = 1;
    HSTR_LOG_DOM log_domain = device / partitions + 1;
    LOG_HSTR_RESULT(hStreams_Alloc1DEx(ptr, num*sizeof(cx_double), NULL, 1, &log_domain));
    LOG_HSTR_RESULT(hStreams_app_xfer_memory(get_stream_id(H2D), ptr, ptr, num*sizeof(cx_double), HSTR_SRC_TO_SINK, NULL));
#   else
    LOG_CUDA_ERROR(cudaSetDevice(device));
    cx_double * ptr;
    LOG_CUDA_ERROR(cudaMallocHost((void**)&ptr, num * sizeof(cx_double)));
    cx_double * device_ptr;
    LOG_CUDA_ERROR(cudaMalloc((void**)&device_ptr, num * sizeof(cx_double)));
    host_to_device[ptr] = device_ptr;
    LOG_CUDA_ERROR(cudaEventCreate(&signals[ptr]));
#   endif
    return ptr;
  }
  void dealloc(cx_double * ptr) {
#   ifdef HS_HSTREAMS
    LOG_HSTR_RESULT(hStreams_DeAlloc(ptr));
    free(ptr);
#   else
    LOG_CUDA_ERROR(cudaSetDevice(device));
    LOG_CUDA_ERROR(cudaFree(host_to_device[ptr]));
    LOG_CUDA_ERROR(cudaEventDestroy(signals[ptr]));
    LOG_CUDA_ERROR(cudaFreeHost(ptr));
#   endif
  }

  static void init(int * num_devices, size_t * total_mem_per_device) {
    int count;
    size_t total; // bytes
#   ifdef HS_HSTREAMS
    LOG_HSTR_RESULT(hStreams_Cfg_SetLogLevel(HSTR_LOG_LEVEL_FATAL_ERROR));
//    LOG_HSTR_RESULT(hStreams_Cfg_SetLogLevel(HSTR_LOG_LEVEL_DEBUG1));
//    LOG_HSTR_RESULT(hStreams_Cfg_SetLogInfoType(-1));
    size_t partitions = 1;
    LOG_HSTR_RESULT(hStreams_app_init(partitions/*partitions*/, 3/*parallel streams*/));
    bool homog = false;
    uint32_t count_uint, nactive, nthr, coreMaxMHz;
    HSTR_ISA_TYPE isa;
    uint64_t suppMemTypes;
    HSTR_CPU_MASK maxMask, avoidMask;
    LOG_HSTR_RESULT(hStreams_GetNumPhysDomains(&count_uint, &nactive, &homog));
    count = count_uint * partitions;
    assert(homog);
    uint64_t mem_info[HSTR_MEM_TYPE_SIZE];
    LOG_HSTR_RESULT(hStreams_GetPhysDomainDetails(0, &nthr, &isa, &coreMaxMHz, maxMask, avoidMask, &suppMemTypes, mem_info));
    // no idea if this is correct once we get offload-over-fabric and other crazy shit.
    total = mem_info[HSTR_MEM_TYPE_HBW] ? mem_info[HSTR_MEM_TYPE_HBW] : mem_info[HSTR_MEM_TYPE_NORMAL];
    total /= partitions;
#   else
    LOG_CUDA_ERROR(cudaGetDeviceCount(&count));
    LOG_CUDA_ERROR(cudaMemGetInfo(nullptr, &total));
    // lets just assume all cuda devices have the same memory.
    // Later: Do check this!
    for (int i = 0; i < count; i++) {
      size_t total_dev;
      LOG_CUDA_ERROR(cudaSetDevice(i));
      LOG_CUDA_ERROR(cudaMemGetInfo(nullptr, &total_dev));
      assert(total_dev == total);
    }
#   endif
    if (num_devices) *num_devices = count;
    if (total_mem_per_device) *total_mem_per_device = total;
  }

  static void release() {
#   ifdef HS_HSTREAMS
    LOG_HSTR_RESULT(hStreams_app_fini());
#   endif
  }
};

struct hs_buffer {
  cx_double * const A, * const B, * const C;
  cx_double * host_C;
  int host_ldC;
  int M, N, K;
  int kind;

  hs_buffer(cx_double * cA, cx_double * cB, cx_double * cC) : A(cA), B(cB), C(cC) {}

  void unpack_C() {
    hs_timing_guard tg(hs_time_unpack);
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
      int lo = i;
      if (kind == HS_GEMM) lo = 0;
      //double * restrict dst = (double *) &host_C[i * host_ldC + lo];
      //double * restrict src = (double *) &C[i * M + lo];
      //const int len = 2 * (M - lo);
      //#pragma omp simd
      //for (int j = 0; j < len; j++) {
      //  dst[i] += src[i];
      //}
      for (int j = lo; j < M; j++) {
        //host_C[i * host_ldC + j].re += C[i * M + j].re;
        //host_C[i * host_ldC + j].im += C[i * M + j].im;
        host_C[i * host_ldC + j] += C[i * M + j];
      }
    }
  }

  void set_params(int theKind, int theM, int theN, int theK, cx_double * theC, int ldC) {
    kind = theKind;
    M = theM;
    N = theN;
    K = theK;
    host_C = theC;
    host_ldC = ldC;
  }

};

struct hs_buffers {
  hs_accelerator accel;
  std::queue<hs_buffer *> empty, transferring, computing;
  int buffer_num;
  bool transferred_one;

  hs_buffers(int dev, int ndevs, int bs, int num_buffers) : accel(dev, ndevs), buffer_num(num_buffers), transferred_one(false) {
    cx_double * A, * B, * C;
    for (int i = 0; i < buffer_num; i++) {
      A = accel.alloc(bs * bs);
      B = accel.alloc(bs * bs);
      C = accel.alloc(bs * bs);
      empty.push(new hs_buffer(A, B, C));
    }
    assert(buffer_num > 0);
    accel.wait(C);
    accel.thunk(HS_GEMM, bs, bs, bs, A, B, C);
    accel.wait(C);
  }

  ~hs_buffers() {
    assert(transferring.empty());
    assert(computing.empty());
    assert(empty.size() == buffer_num);
    for (int i = 0; i < buffer_num; i++) {
      hs_buffer * buf = empty.front();
      empty.pop();
      accel.dealloc(buf->A);
      accel.dealloc(buf->B);
      accel.dealloc(buf->C);
      delete buf;
    }
  }

  hs_buffer * get_first_empty_buffer(bool host_compute) {
    advance(host_compute ? 0 : 1);
    //assert(!empty.empty());
    if (empty.empty()) {
      return nullptr;
    } else {
      return empty.front();
    }
  }

  void issue_async_computes() {
    int trans_size = transferring.size();
    for (int i = 0; i < trans_size; i++) {
      hs_buffer * buf = transferring.front();
      transferring.pop();
      if (accel.check(buf->A) && accel.check(buf->B)) {
        accel.thunk(buf->kind, buf->M, buf->N, buf->K, buf->A, buf->B, buf->C);
        computing.push(buf);
        trans_size -= 1;
        transferred_one = true;
      } else {
        transferring.push(buf);
      }
    }
  }

  void send_compute() {
    if (! transferred_one && ! transferring.empty()) {
      hs_buffer * done = transferring.front();
      accel.thunk(done->kind, done->M, done->N, done->K, done->A, done->B, done->C);
      transferring.pop();
      computing.push(done);
    }
    transferred_one = false;
  }

  void send_all_computes() {
    while (! transferring.empty()) {
      send_compute();
    }
  }

  void advance(int threshold) {
    send_compute();
    if ((! computing.empty()) && (empty.size() < threshold)) {
      hs_buffer * done = computing.front();
      accel.wait(done->C);
      done->unpack_C();
      computing.pop();
      empty.push(done);
    } else { // threshold is 0
      int comp_size = computing.size();
      for (int i = 0; i < comp_size; i++) {
        hs_buffer * done = computing.front();
        computing.pop();
        if (accel.check(done->C)) {
          done->unpack_C();
          empty.push(done);
          comp_size -= 1;
        } else {
          computing.push(done);
        }
      }
    }
  }

  void enqueue_first_empty_buffer(hs_buffer * slot) {
    assert(!empty.empty());
    hs_buffer * first = empty.front();
    assert(first == slot);
    empty.pop();
    transferring.push(first);
  }

  const cx_double * pack_A(hs_buffer * at, const cx_double * buf, int ld, int num_strides, int num_in_stride, bool try_avoid_buffer) {
    return accel.transfer_strided(at->A, buf, ld, num_strides, num_in_stride, try_avoid_buffer);
  }

  const cx_double * pack_B(hs_buffer * at, const cx_double * buf, int ld, int num_strides, int num_in_stride, bool try_avoid_buffer) {
    return accel.transfer_strided(at->B, buf, ld, num_strides, num_in_stride, try_avoid_buffer);
  }
};

void hs_kernel_block(int &split_blocks, int myBS, int kind, int i0, int j0, int k0, int M, int N, int K, const cx_double * A, int ldA, const cx_double * B, int ldB, double beta, cx_double * C, int ldC);
static inline int max(int a, int b) { return a > b ? a : b; }

struct hs_manager {
  int the_device;
  int num_devices;
  int the_BS;
  bool host_compute;
  int num_buffers;
  int split_blocks;
  bool try_avoid_buffer;
  bool is_inside_host_recursion;
  std::vector< std::unique_ptr<hs_buffers> > buffers;

  hs_manager() : num_devices(0), the_device(0), the_BS(1000), is_inside_host_recursion(false) {}

  int optional_env_int(char * str, int orig) {
    char * env = getenv(str);
    if (env) {
      printf("%s %d -> %d\n", str, orig, atoi(env));
      orig = atoi(env);
    }
    return orig;
  }

  void init(size_t mem, int ND) {
    assert(num_devices == 0);
    num_buffers = optional_env_int("HSNB", 5);
    num_devices = optional_env_int("HSND", ND);
    host_compute = optional_env_int("HSHC", true);
    split_blocks = optional_env_int("HSSPLIT", 0);
    int mul = optional_env_int("HSMUL", 1);
    int BS = sqrt(0.6 * mem / 3 / num_buffers / mul/ sizeof(double) / 2);
    int factor = optional_env_int("HSNDSTRIDE", 1);
    the_BS = optional_env_int("HSBS", BS);
    try_avoid_buffer = optional_env_int("HSTAVB", 1);
    // HSHR=0 disables host recursion
    is_inside_host_recursion = ! optional_env_int("HSHR", 1);
    for (int j = 0; j < mul; j++) {
      for (int i = 0; i < num_devices; i++) {
        buffers.emplace_back(new hs_buffers(factor * i, ND, BS, num_buffers));
      }
    }
    printf("[DEVICES] %d devices\n", num_devices );
  }

  void release() {
    buffers.clear();
    num_devices = 0;
    the_device = 0;
  }

  void begin_calculation() {
    assert(the_device == 0);
#   ifdef HS_CUDA
    cudaSetDevice(the_device);
#   endif
  }

  hs_buffers * get_fresh_buffer() {
    return buffers[the_device++ % buffers.size()].get();
  }

  void end_calculation() {
    for (int i = 0; i < buffers.size() * num_buffers * 2; i++) {
      hs_buffers * buf = get_fresh_buffer();
      buf->transferred_one = false;
      buf->send_all_computes();
    }
    for (int i = 0; i < buffers.size() * num_buffers * 2; i++) {
      hs_buffers * buf = get_fresh_buffer();
      buf->transferred_one = false;
      buf->advance(num_buffers);
    }
    the_device = 0;
#   ifdef HS_CUDA
    cudaSetDevice(the_device);
#   endif
  }

  void enqueue_compute(int kind, int theM, int theN, int theK, const cx_double * A, int ldA, const cx_double * B, int ldB, cx_double * C, int ldC) {
    hs_buffers * buf = nullptr;
    hs_buffer * slot = nullptr;
    for (int i = 0; i < buffers.size(); i++) {
      buffers[i]->issue_async_computes();
    }
    for (int i = 0; i < buffers.size(); i++) {
      buf = get_fresh_buffer();
      slot = buf->get_first_empty_buffer(host_compute);
      if (slot) break;
    }
    if (slot) {
      buf->pack_A(slot, A, ldA, theM, theK, try_avoid_buffer);
      if (kind == HS_HERK) {
        buf->pack_B(slot, A, 1, 1, 1, try_avoid_buffer);
      } else {
        buf->pack_B(slot, B, ldB, theN, theK, try_avoid_buffer);
      }
      slot->set_params(kind, theM, theN, theK, C, ldC);
      buf->enqueue_first_empty_buffer(slot);
    } else if (! is_inside_host_recursion) {
      is_inside_host_recursion = true;
      int split_is_zero = 0;
      hs_kernel_block(split_is_zero, max(theM, max(theN, theK)) / 2, kind, 0, 0, 0, theM, theN, theK, A, ldA, B, ldB, 1, C, ldC);
      is_inside_host_recursion = false;
    } else {
      hs_timing_guard tg(hs_host_compute);
      //cx_double c0 = {0}, c1 = {0};
      //c1.re = 1;
      cx_double c0(0, 0), c1(1, 0);
      double d1 = 1.0;
      if (kind == HS_GEMM) {
        zgemm_("C", "N", &theM, &theN, &theK, &c1, A, &ldA, B, &ldB, &c1, C, &ldC);
      } else if (kind == HS_HERK) {
        zherk_("L", "C", &theN, &theK, &d1, A, &ldA, &d1, C, &ldC);
      } else if (kind == HS_HER2K) {
        zher2k_("L", "C", &theN, &theK, &c1, A, &ldA, B, &ldB, &d1, C, &ldC);
      } else if (kind == HS_HERKX) {
        zgemmt_("L", "C", "N", &theN, &theK, &c1, A, &ldA, B, &ldB, &c1, C, &ldC);
        for (int i = 0; i < theN; i++) { // prob not worth parallelizing.
          cx_double val(C[i * ldC + i].real(), 0);
          C[i * ldC + i] = val;
        }
      }
    }
  }
} hs_g_manager;

void hs_init(void * A, void * B, uint64_t szAB, void * H, void * S, uint64_t szHS) {
  int nd;
  size_t mem;
  hs_accelerator::init(&nd, &mem);
  hs_g_manager.init(mem, nd);
}

void hs_kernel_block(int &split_blocks, int myBS, int kind, int i0, int j0, int k0, int M, int N, int K, const cx_double * A, int ldA, const cx_double * B, int ldB, double beta, cx_double * C, int ldC) {
  //cx_double c0 = {0}, c1 = {0};
  //c1.re = 1;
  cx_double c0(0, 0), c1(1, 0);
  double d1 = 1.0;
  for (int i = i0; i < N; i += myBS) { // column of C
    int j_lo = i; // on the floor
    if (kind == HS_GEMM) j_lo = 0;
    if (j_lo < j0) j_lo = j0;
    for (int j = j_lo; j < M; j += myBS) { // row of C
      for (int k = k0; k < K; k += myBS) {
        int theN = myBS;
        if (i + myBS > N) theN = N - i;
        int theM = myBS;
        if (j + myBS > M) theM = M - j;
        int theK = myBS;
        if (k + myBS > K) theK = K - k;
        if ((beta == 0.0) && (k == 0)) {
          hs_timing_guard tg(hs_time_zero);
          // zero the C-block out.
          #pragma omp parallel for
          for (int ii = 0; ii < theN; ii++) {
            memset(&C[(i + ii) * ldC + j], 0, theM * sizeof(cx_double));
          }
        }
        if (split_blocks) {
          split_blocks -= 1;
          hs_kernel_block(split_blocks, (myBS + 1)/2, kind, i, j, k, j+theM, i+theN, k+theK, A, ldA, B, ldB, 1, C, ldC);
          continue;
        }
        if (i == j) {
          hs_g_manager.enqueue_compute(kind, theM, theN, theK, &A[j * ldA + k], ldA, &B[i * ldB + k], ldB, &C[i * ldC + j], ldC);
        } else {
          hs_g_manager.enqueue_compute(HS_GEMM, theM, theN, theK, &A[j * ldA + k], ldA, &B[i * ldB + k], ldB, &C[i * ldC + j], ldC);
          if (kind == HS_HER2K) {
            hs_g_manager.enqueue_compute(HS_GEMM, theM, theN, theK, &B[j * ldB + k], ldB, &A[i * ldA + k], ldA, &C[i * ldC + j], ldC);
          }
        }
      }
    }
  }
}

int hs_pick_bs(int BS, int M, int N, int K) {
  // requirements:
  // * do not exceed the_BS.
  // * have it as evenly as possible split the smallest dimension
  // * (maybe) have it such that the remainder is minimized in dimensions proportional to their length
  // * (maybe) have at least two blocks in each dimension
  // * divide by 64/32
  // * do not divide by 2048/4096
  int dims[3] = {M, N, K};
  int bss[3] = {0};
  for (int i = 0; i < 3; i++) {
    int dim_min = dims[i];
    int nb_min = (dim_min + BS - 1) / BS; // ceil(dim_min / BS)
    if (nb_min == 1) nb_min += 1;
    int bs_min = (dim_min + nb_min - 1) / nb_min; // ceil(dim_min / nb_min)
    assert(bs_min <= BS);
    const int ELEM_MUL = 64;
    while (bs_min % ELEM_MUL != 0) bs_min += 1;
    if (bs_min > BS) bs_min -= ELEM_MUL;
    const int ELEM_NOT_MUL = 256;
    if (bs_min % ELEM_NOT_MUL == 0) bs_min += ELEM_MUL;
    if (bs_min > BS) bs_min -= 2 * ELEM_MUL;
    bss[i] = bs_min;
  }
  int score[3] = {1, 1, 1};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      score[i] *= dims[j] % bss[i];
    }
  }
  int max_score_idx = 0;
  for (int i = 1; i < 3; i++) {
    if (score[i] > score[max_score_idx]) max_score_idx = i;
  }
  int bs_min = bss[max_score_idx];
  printf("PICK %d -> %d %d %d -> %d\n", BS, M, N, K, bs_min);
  return bs_min;
}

void hs_kernel(int kind, int M, int N, int K, const cx_double * A, int ldA, const cx_double * B, int ldB, double beta, cx_double * C, int ldC) {
  hs_timing_guard tg(hs_time_total);
  //#pragma omp parallel num_threads(2)
  //{
  //  if (omp_get_thread_num() == 0) {
  //    omp_set_num_threads(4);
  //    omp_set_nested(1);
  //    #pragma omp parallel
  //    {
  //      kmp_affinity_mask_t mask;
  //      kmp_create_affinity_mask(&mask);
  //      kmp_set_affinity_mask_proc(omp_get_thread_num(), &mask);
  //      if (kmp_set_affinity(&mask) != 0)
  //       printf("Could not set affinity on master %d.\n", omp_get_thread_num());
  //    }
      int myBS = hs_pick_bs(hs_g_manager.the_BS, M, N, K);
      hs_g_manager.begin_calculation();
      int split = hs_g_manager.split_blocks;
      hs_kernel_block(split, myBS, kind, 0, 0, 0, M, N, K, A, ldA, B, ldB, beta, C, ldC);
      hs_g_manager.end_calculation();
  //  } else if (omp_get_thread_num() == 1) {
  //    omp_set_num_threads(12);
  //    omp_set_nested(1);
  //    #pragma omp parallel
  //    {
  //      kmp_affinity_mask_t mask;
  //      kmp_create_affinity_mask(&mask);
  //      kmp_set_affinity_mask_proc(4 + omp_get_thread_num(), &mask);
  //      if (kmp_set_affinity(&mask) != 0)
  //       printf("Could not set affinity on worker %d.\n", omp_get_thread_num());
  //    }
  //    // do the work
  //  } else { abort(); }
  //}
}

void hs_zher2k(const char * uplo, const char * trans, 
               const int * N, const int * K, 
               const cx_double * alpha, const cx_double * A, const int * ldA, 
                                        const cx_double * B, const int * ldB, 
               const double * beta, cx_double * C, const int * ldC) {
  //if (*uplo != 'L' || *trans != 'C' || *alpha != cx_double(1,0) || *beta != 0.0 || *N == 0 || *K == 0) {
  if (*uplo != 'L' || *trans != 'C' || fabs(alpha->real()-1.0) > 1e-15 || fabs(alpha->imag()) > 1e-15 ||
                                       fabs(*beta) > 1e-15 || *N == 0 || *K == 0) {
    printf("Warning: Running on host.\n");
    zher2k_(uplo, trans, N, K, alpha, A, ldA, B, ldB, beta, C, ldC);
    return;
  }
//  bool old_try = hs_g_manager.try_avoid_buffer;
  //hs_g_manager.try_avoid_buffer = true;
  hs_kernel(HS_HER2K, *N, *N, *K, A, *ldA, B, *ldB, 0, C, *ldC);
  //hs_g_manager.try_avoid_buffer = old_try;
}

void hs_zherk(const char * uplo, const char * trans, 
              const int * N, const int * K, 
              const double * alpha, const cx_double * A, const int * ldA, 
              const double * beta,        cx_double * C, const int * ldC) {
  //if (*uplo != 'L' || *trans != 'C' || fabs(*alpha != 1.0 || (*beta != 0.0 && *beta != 1.0) || *N == 0 || *K == 0) {
  if (*uplo != 'L' || *trans != 'C' || fabs(*alpha-1.0) > 1e-15 || (fabs(*beta) > 1e-15 && fabs(*beta-1.0) > 1e-15) || *N == 0 || *K == 0) {
    printf("Warning: Running on host.\n");
    zherk_(uplo, trans, N, K, alpha, A, ldA, beta, C, ldC);
    return;
  }
  hs_kernel(HS_HERK, *N, *N, *K, A, *ldA, A, *ldA, *beta, C, *ldC);
}

void hs_zgemm(const char * transA, const char * transB, 
              const int * M, const int * N, const int * K, 
              const cx_double * alpha, const cx_double * A, const int * ldA, 
                                       const cx_double * B, const int * ldB, 
              const cx_double * beta,        cx_double * C, const int * ldC) {
  //if (*transA != 'C' || *transB != 'N' || *alpha != cx_double(1,0) || (*beta != cx_double(1,0) && *beta != cx_double(0,0)) || *M == 0 || *N == 0 || *K == 0) {
  if (*transA != 'C' || *transB != 'N' || fabs(alpha->real()-1.0) > 1e-15 ||  fabs(alpha->imag()) > 1e-15 || 
                                         (fabs(beta->real()-1.0) > 1e-15 && fabs(beta->real()) > 1e-15) || 
                                          fabs(beta->imag()) > 1e-15 || *N == 0 || *K == 0) {
    printf("Warning: Running on host.\n");
    zgemm_(transA, transB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC);
    return;
  }
  hs_kernel(HS_GEMM, *M, *N, *K, A, *ldA, B, *ldB, beta->real(), C, *ldC);
}

void hs_zherkx(const char * uplo, const char * trans, 
               const int * N, const int * K, 
               const cx_double * alpha, const cx_double * A, const int * ldA, 
                                        const cx_double * B, const int * ldB, 
               const double * beta,           cx_double * C, const int * ldC) {
  //if (*uplo != 'L' || *trans != 'C' || *alpha != cx_double(1,0) || (*beta != 1.0d && *beta != 0.0d) || *N == 0 || *K == 0) {
  if (*uplo != 'L' || *trans != 'C' || fabs(alpha->real()-1.0) > 1e-15 || fabs(alpha->imag()) > 1e-15 ||
                                      (fabs(*beta) > 1e-15 && fabs(*beta-1.0) > 1e-15) || *N == 0 || *K == 0) {
    printf("Warning: Running on host.\n");
    cx_double cx_beta(*beta,0);
    zgemmt_(uplo, trans, *trans == 'C' ? "N" : "C", N, K, 
            alpha, A, ldA, B, ldB, &cx_beta, C, ldC);
    return;
  }
  hs_kernel(HS_HERKX, *N, *N, *K, A, *ldA, B, *ldB, *beta, C, *ldC);
}

void hs_fini() {
  hs_g_manager.release();
  printf("Timings---------\n");
  hs_wait_xfer.output("wait xfer");
  hs_wait_compute.output("wait comp");
  hs_wait_issue_xfer.output("wais xfer");
  hs_wait_issue_compute.output("wais comp");
  hs_wait_issue_event.output("wais evnt");
  hs_wait_issue_output.output("wais outp");
  hs_time_unpack.output("time unpa");
  hs_time_pack.output("time pack");
  hs_time_zero.output("time zero");
  hs_check_compute.output("chek comp");
  hs_host_compute.output("host comp");
  hs_time_total.output("time totl");
  hs_accelerator::release();
}
#endif

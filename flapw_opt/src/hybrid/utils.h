/*
 * test_utils.h
 *
 *  Created on: Sep 26, 2012
 *      Author: ddavidov
 */


#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuda.h>

#ifndef UTILS_H_
#define UTILS_H_

#define LOG_CUDA_ERROR(func)                                                 \
{                                                                         \
  cublasStatus_t hret = CUBLAS_STATUS_SUCCESS;                               \
  hret = func;                                                          \
  if (hret != CUBLAS_STATUS_SUCCESS) {                                    \
     printf("%s returned %s. (%s:%d)\n", #func, cublasGetErrorString(hret), __FILE__, __LINE__); \
     abort();                                                           \
  }                                                                     \
}

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

#define LOG_CUBLAS_STATUS(func)                                                 \
{                                                                         \
  cublasStatus_t hret = CUBLAS_STATUS_SUCCESS;                               \
  hret = func;                                                          \
  if (hret != CUBLAS_STATUS_SUCCESS) {                                    \
    printf("%s returned %s. (%s:%d)\n", #func, cublasGetErrorString(hret), __FILE__, __LINE__); \
    abort();                                                           \
  }                                                                     \
}

#endif /* UTILS_H_ */

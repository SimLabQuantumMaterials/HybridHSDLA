/*
 * utils.c
 *
 *  Created on: June 10, 2018
 *      Author: Davor Davidovic, Rudjer Boskovic Institute, Croatia
 */


#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuda.h>

#include "utils.h"

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

void LOG_CUDA_STATUS(cudaError_t hret, char *func )
{
  if (hret != cudaSuccess){
     printf("%s returned %s. (%s:%d)\n", func, cudaGetErrorString(hret), __FILE__, __LINE__);
     abort();
  }
}


void LOG_CUBLAS_STATUS(cublasStatus_t hret, char *func)
{
  if (hret != CUBLAS_STATUS_SUCCESS){
    printf("%s returned %s. (%s:%d)\n", func, cublasGetErrorString(hret), __FILE__, __LINE__);
    abort();
  }
}

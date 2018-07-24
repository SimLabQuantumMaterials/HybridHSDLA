/*
 * utils.h
 *
 *  Created on: June 10, 2018
 *      Author: Davor Davidovic, Rudjer Boskovic Institute, Croatia
 */


#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuda.h>

const char * cublasGetErrorString(cublasStatus_t);

void LOG_CUDA_STATUS(cudaError_t, char * );

void LOG_CUBLAS_STATUS(cublasStatus_t, char *);

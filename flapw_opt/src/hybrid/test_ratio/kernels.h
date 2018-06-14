/*
 *  	Definition of the function to compute ratios
 */

#ifndef KERNELS_H_
#define KERNELS_H_

#include <complex>
#include <cublasXt.h>

/* ZHERK */
void ratio_zherk( cublasXtHandle_t, int, int, int, int, std::complex<double>*, int, std::complex<double>*, int, double* );

/* ZHERKX */
void ratio_zherkx( cublasXtHandle_t, int, int, int, int, std::complex<double>*, int, std::complex<double>*, int, std::complex<double>*, int, double* );

/* ZHER2K */
void ratio_zher2k( cublasXtHandle_t, int, int, int, int, std::complex<double>*, int, std::complex<double>*, int, std::complex<double>*, int, double* );

#endif

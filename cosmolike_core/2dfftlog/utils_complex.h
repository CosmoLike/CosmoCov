/******************************************************************************
adapted from 2DFFTLog
https://github.com/xfangcosmo/2DFFTLog
by Xiao Fang
******************************************************************************/

#include <complex.h>
#include <fftw3.h>

void g_l_cfastcov(double l, double nu, double *eta, double complex *gl, long N);

void g_l_smooth_cfastcov(double l, double nu, double *eta, double complex *gl, long N, double smooth_dlnr, double alpha_pow);

void c_window_2d(double complex *out, double c_window_width, long halfN1, long halfN2);

// void resample_fourier_gauss(double *k, double *fk, config *config);

double complex gamma_lanczos_cfastcov(double complex z);
double complex lngamma_lanczos_cfastcov(double complex z);
/******************************************************************************
adapted from 2DFFTLog
https://github.com/xfangcosmo/2DFFTLog
by Xiao Fang
******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "utils.h"

double complex gamma_lanczos_cfastcov(double complex z) {
/* Lanczos coefficients for g = 7 */
	static double p[] = {
		0.99999999999980993227684700473478,
		676.520368121885098567009190444019,
		-1259.13921672240287047156078755283,
		771.3234287776530788486528258894,
		-176.61502916214059906584551354,
		12.507343278686904814458936853,
		-0.13857109526572011689554707,
		9.984369578019570859563e-6,
		1.50563273514931155834e-7};

	if(creal(z) < 0.5) {return M_PI / (csin(M_PI*z)*gamma_lanczos_cfastcov(1. - z));}
	z -= 1;
	double complex x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double complex t = z + 7.5;
	return sqrt(2*M_PI) * cpow(t, z+0.5) * cexp(-t) * x;
}


double complex lngamma_lanczos_cfastcov(double complex z) {
/* Lanczos coefficients for g = 7 */
	static double p[] = {
		0.99999999999980993227684700473478,
		676.520368121885098567009190444019,
		-1259.13921672240287047156078755283,
		771.3234287776530788486528258894,
		-176.61502916214059906584551354,
		12.507343278686904814458936853,
		-0.13857109526572011689554707,
		9.984369578019570859563e-6,
		1.50563273514931155834e-7};

	if(creal(z) < 0.5) {return log(M_PI) -clog(csin(M_PI*z)) - lngamma_lanczos_cfastcov(1. - z);}
	z -= 1;
	double complex x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	double complex t = z + 7.5;
	return log(2*M_PI) /2.  + (z+0.5)*clog(t) -t + clog(x);
}

void g_l_cfastcov(double l, double nu, double *eta, double complex *gl, long N) {
/* z = nu + I*eta
Calculate g_l = zln2 + lngamma( (l+nu)/2 + I*eta/2 ) - lngamma( (3+l-nu)/2 - I*eta/2 ) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl[i] = cexp(z*log(2.) + lngamma_lanczos_cfastcov((l+z)/2.) - lngamma_lanczos_cfastcov((3.+l-z)/2.));
	}
}

void g_l_smooth_cfastcov(double l, double nu, double *eta, double complex *gl, long N, double smooth_dlnr, double alpha_pow) {
/* z = nu + I*eta
Calculate g_l_smooth = g_l * exp((2.-z)*smooth_dlnr) / (2.-z) */
	long i;
	double complex z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl[i] = cexp(z*log(2.) + lngamma_lanczos_cfastcov((l+z)/2.) - lngamma_lanczos_cfastcov((3.+l-z)/2.));
		gl[i] *= (cexp((alpha_pow-z)*smooth_dlnr)-1.)/(alpha_pow-z);
	}
}

void c_window_2d(double complex *out, double c_window_width, long halfN1, long halfN2) {
	// 'out' is N1*(halfN2+1) complex array
	long Ncut1, Ncut2;
	long N1;
	N1 = 2*halfN1;
	Ncut1 = (long)(halfN1 * c_window_width);
	Ncut2 = (long)(halfN2 * c_window_width);
	long i,j;
	double W;
	for(j=0; j<=Ncut2; j++) { // window for right-side
		W = (double)(j)/Ncut2 - 1./(2.*M_PI) * sin(2.*j*M_PI/Ncut2);
		for(i=0; i<N1; i++){
			out[i*(halfN2+1)+ halfN2-j] *= W;
		}
	}
	for(i=0; i<=Ncut1; i++){ // window for center part (equivalent for up-down edges when rearanged)
		W = (double)(i)/Ncut1 - 1./(2.*M_PI) * sin(2.*i*M_PI/Ncut1);
		for(j=0; j<=halfN2; j++) {
			out[(halfN1-i)*(halfN2+1) + j] *= W;
			out[(halfN1+1+i)*(halfN2+1) + j] *= W;
		}
	}
}


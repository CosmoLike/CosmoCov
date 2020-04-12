/******************************************************************************
adapted from 2DFFTLog
https://github.com/xfangcosmo/2DFFTLog
by Xiao Fang
******************************************************************************/

#include <stdlib.h>
#include <math.h>

#include "utils.h"


void mk_diag_g_to_ng(double *in, long N, double dlnk, double **out) {
	long i, j;
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			out[i][j] = 0.;
		}
		out[i][i] = in[i] / dlnk;
	}
}

void extrap_log_linear_cfastcov(double *fk, int N_origin, int N_extra, double *large_fk) {
	double dln_left, dln_right;
	int i;

	dln_left = log(fk[1]/fk[0]);
	// printf("fk[0],fk[1]: %.15e,%.15e,%.15e,%.15e,%.15e\n", fk[0],fk[1],fk[2],fk[3],fk[4]);
	if(fk[0]<=0.) {
		for(i=0; i<N_extra; i++) {
			large_fk[i] = 0.;
		}
	}
	else{
		for(i=0; i<N_extra; i++) {
			large_fk[i] = exp(log(fk[0]) + (i - N_extra) * dln_left);
		}
	}

	for(i=N_extra; i< N_extra+N_origin; i++) {
		large_fk[i] = fk[i - N_extra];
	}

	dln_right = log(fk[N_origin-1]/fk[N_origin-2]);
	if(fk[N_origin-1]<=0.) {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = 0.;
		}
	}
	else {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = exp(log(fk[N_origin-1]) + (i - N_extra - N_origin +1) * dln_right);
		}
	}
}

void extrap_log_bilinear(double **fk, int N_origin, int N_extra, double **large_fk) {
	int i,j;
	double dln_left, dln_right;
	for(i=N_extra; i<N_origin+N_extra; i++) {
		dln_left = log(fk[i-N_extra][1]/fk[i-N_extra][0]);
		dln_right = log(fk[i-N_extra][N_origin-1]/fk[i-N_extra][N_origin-2]);
		for(j=0; j<N_extra; j++) {
			large_fk[i][j] = exp(log(fk[i-N_extra][0]) + (j - N_extra) * dln_left);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(j=N_extra; j< N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][j - N_extra];
		}
		for(j=N_extra+N_origin; j< 2*N_extra+N_origin; j++) {
			large_fk[i][j] = exp(log(fk[i-N_extra][N_origin-1]) + (j - N_extra - N_origin +1) * dln_right);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
	double dln_up, dln_down;
	for(j=0; j<N_origin+2*N_extra; j++) {
		dln_up = log(large_fk[N_extra+1][j]/large_fk[N_extra][j]);
		dln_down = log(large_fk[N_extra+N_origin-1][j]/large_fk[N_extra+N_origin-2][j]);
		for(i=0; i<N_extra; i++) {
			large_fk[i][j] = exp(log(large_fk[N_extra][j]) + (i - N_extra) * dln_up);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i][j] = exp(log(large_fk[N_extra+N_origin-1][j]) + (i - N_extra - N_origin +1) * dln_down);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
}


void extrap_bilinear(double **fk, int N_origin, int N_extra, double **large_fk) {
	int i,j;
	double dleft, dright;
	for(i=N_extra; i<N_origin+N_extra; i++) {
		dleft = fk[i-N_extra][1]-fk[i-N_extra][0];
		dright = fk[i-N_extra][N_origin-1]-fk[i-N_extra][N_origin-2];
		for(j=0; j<N_extra; j++) {
			large_fk[i][j] = fk[i-N_extra][0] + (j - N_extra) * dleft;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(j=N_extra; j< N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][j - N_extra];
		}
		for(j=N_extra+N_origin; j< 2*N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][N_origin-1] + (j - N_extra - N_origin +1) * dright;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
	double dup, ddown;
	for(j=0; j<N_origin+2*N_extra; j++) {
		dup = large_fk[N_extra+1][j]-large_fk[N_extra][j];
		ddown = large_fk[N_extra+N_origin-1][j]-large_fk[N_extra+N_origin-2][j];
		for(i=0; i<N_extra; i++) {
			large_fk[i][j] = large_fk[N_extra][j] + (i - N_extra) * dup;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i][j] = large_fk[N_extra+N_origin-1][j] + (i - N_extra - N_origin +1) * ddown;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
}

void extrap_2dzeros(double **fk, int N_origin, int N_extra, double **large_fk) {
	int i,j;
	for(i=N_extra; i<N_origin+N_extra; i++) {
		for(j=0; j<N_extra; j++) {
			large_fk[i][j] = 0.;
		}
		for(j=N_extra; j< N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][j - N_extra];
			// if(isnan(large_fk[i][j])) {printf("large_fk nan at: %d,%d\n", i,j);large_fk[i][j]=0.;}
		}
		for(j=N_extra+N_origin; j< 2*N_extra+N_origin; j++) {
			large_fk[i][j] = 0.;
		}
	}
	for(j=0; j<N_origin+2*N_extra; j++) {
		for(i=0; i<N_extra; i++) {
			large_fk[i][j] = 0.;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i][j] = 0.;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
}

// void resample_fourier_gauss(double *k, double *fk, config *config, double *k_sample, double *fk_sample) {
// 	long i;
// 	double dlnk = log(k[sizeof(k)-1]/k[0]) / (config->Nk_sample-1.);
// 	for(i=0; i<config->Nk_sample; i++) {
// 		k_sample[i] = k[0] * exp(i*dlnk);
// 		fk_sample[i] = 
// 	}
// }
/******************************************************************************
adapted from 2DFFTLog
https://github.com/xfangcosmo/2DFFTLog
by Xiao Fang
******************************************************************************/

typedef struct config {
	int sym_Flag;
	double l1, l2;
	double nu1, nu2;
	double c_window_width;
	long Nk_sample;
} config;
void twobessel(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double *r1, double *r2, double **result);
void twobessel_binave(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double smooth_dlnr, double alpha_pow, double *r1, double *r2, double **result);
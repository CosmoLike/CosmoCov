/******************************************************************************
adapted from 2DFFTLog
https://github.com/xfangcosmo/2DFFTLog
by Xiao Fang
******************************************************************************/

typedef struct preconfig {
	double pre_k1min, pre_k1max;
	double pre_k2min, pre_k2max;
} preconfig;

void mk_diag_g_to_ng(double *in, long N, double dlnk, double **out);

void extrap_log_linear_cfastcov(double *fk, int N_origin, int N_extra, double *large_fk);

void extrap_log_bilinear(double **fk, int N_origin, int N_extra, double **large_fk);

void extrap_bilinear(double **fk, int N_origin, int N_extra, double **large_fk);

void extrap_2dzeros(double **fk, int N_origin, int N_extra, double **large_fk);
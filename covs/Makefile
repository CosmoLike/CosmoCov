cfastcov_dir := ../cosmolike_core/2dfftlog/
cfastcov := $(cfastcov_dir)twobessel.c $(cfastcov_dir)utils.c $(cfastcov_dir)utils_complex.c
opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99

cfftlog_dir := ../cosmolike_core/cfftlog/
cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c

covs:
	make cov_fft
	make cov_curved

cov_fft:
	gcc -o cov_flat_fft compute_covariances_real_flat_fft.c $(cfastcov) $(opt_home)

cov_curved:
	gcc -o cov compute_covariances_real.c $(cfftlog) $(opt_home)

response:
	gcc -o cov_response compute_covariances_real_WL_response.c  $(opt_home) -lpython2.7


opt_nersc := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I ${GSL_DIR}/include  -L ${GSL_DIR}/lib -I ${FFTW_INC} -L ${FFTW_DIR} -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99

nersc_covs:
	make nersc_cov_fft
	make nersc_cov_curved

nersc_cov_fft:
	gcc -o cov_flat_fft compute_covariances_real_flat_fft.c $(cfastcov) $(opt_nersc)

nersc_cov_curved:
	gcc -o cov compute_covariances_real.c $(cfftlog) $(opt_nersc)

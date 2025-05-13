#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

// Standard libraries for C

// Parameters
const double a_init = 1e-4; // Initial scale-factor for growth ODE integration

double _w(double a, double w0, double wa) {
    return w0 + (1.0 - a) * wa;
}

double _X_w(double a, double w0, double wa) {
    return pow(a, -3.0 * (1.0 + w0 + wa)) * exp(-3.0 * wa * (1.0 - a));
}

double _Hubble2(double a, double omega0m, double omega0w, double omega0k, double w0, double wa) {
    // Get Om_c, Om_b, and Om_nu using the CAMB_results object
    double Om_m = omega0m;
    double Om_w = omega0w;
    double Om = 1-omega0k;
    double H2 = Om_m * pow(a, -3.0) + Om_w * _X_w(a, w0, wa) + (1.0 - Om) * pow(a, -2.0);
    return H2;
}
double _Omega_m(double a, double omega0m, double omega0w, double omega0k, double w0, double wa) {
    return omega0m * pow(a, -3.0) / _Hubble2(a, omega0m, omega0w, omega0k, w0, wa);
}


double _AH(double a, double omega0m, double omega0w, double w0, double wa) {
    // Get Om_c, Om_b, and Om_nu using the CAMB_results object
    double Om_m = omega0m;
    double Om_w = omega0w;
    double AH = -0.5 * (Om_m * pow(a, -3.0) + (1.0 + 3.0 * _w(a, w0, wa)) * Om_w * _X_w(a, w0, wa));
    return AH;
}
typedef struct {
    double omega0m;
    double omega0w;
    double w0;
    double wa;
    double omega0k;
} Params;




int meadgrowth(double a, const double y[], double f[], void *params) {
    Params *p = (Params *)params;
    f[0] = y[1];
    double  fv = -(2.+_AH(a, p->omega0m, p->omega0w, p->w0, p->wa)/_Hubble2(a, p->omega0m, p->omega0w, p->omega0k, p->w0, p->wa))*y[1]/a;
    double fd = 1.5*_Omega_m(a, p->omega0m, p->omega0w, p->omega0k, p->w0, p->wa)*y[0]/pow(a,2);
    f[1] = fv+fd;
    return GSL_SUCCESS;
}

int main() {
    Params params = {0.3, 0.7, -1, 0., 0};
    gsl_odeiv2_system sys = {meadgrowth, NULL, 2, &params};

    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double a_init=1E-4, a = 1;
//omega0m, double omega0w, double omega0k, double w0, double wa
    double f = 1.-_Omega_m(a_init, 3.120461e-01,6.878622e-01, 0.0, -1, 0);

    double y[2] = {pow(a_init,(1.-3.*f/5.)), (1.-3.*f/5.)*pow(a_init,(-3.*f/5.))};
    int status = gsl_odeiv2_driver_apply(driver, &a_init, a, y);
    printf("%le \n", a_init);

    if (status != GSL_SUCCESS) {
        printf("Error, return value=%d\n", status);
    }

    printf("t=%.5e, y(t)=%.5e, y'(t)=%.5e\n", a, y[0], y[1]);
    gsl_odeiv2_driver_free(driver);
    return 0;
}


#ifndef APPROX1_H
#define APPROX1_H

#ifdef __cplusplus
extern "C" {
#endif

void build_approx1(int n, const double *x, const double *f, double *a, double *work);
double eval_approx1(double t, double xl, double xr, int n, const double *x, const double *a);

#ifdef __cplusplus
}
#endif

#endif
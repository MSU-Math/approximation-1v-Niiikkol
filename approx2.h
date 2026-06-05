#ifndef APPROX2_H
#define APPROX2_H

#ifdef __cplusplus
extern "C" {
#endif

void build_approx2(int n, const double *x, const double *f, double *a, double *work);
double eval_approx2(double t, double xl, double xr, int n, const double *x, const double *a);

#ifdef __cplusplus
}
#endif

#endif 

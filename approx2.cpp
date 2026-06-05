#include "approx2.h"

void build_approx2(int n, const double *x, const double *f, double *a, double *work)
{
    int i;
    double h, bi;
    double hi1;
    double *d = work;
    double *c = work + n; 

    if (n < 2) {
        return;
    }

    for (i = 0; i < n - 1; i++) {
        h = x[i + 1] - x[i];
        d[i] = (f[i + 1] - f[i]) / h;
    }

    if (n == 2) {
        c[0] = 0.0;
    } else {

        c[1] = (d[1] - d[0]) / (x[2] - x[0]);
        c[0] = c[1];

        for (i = 1; i < n - 2; i++) {
            h = x[i + 1] - x[i];
            hi1 = x[i + 2] - x[i + 1];
            c[i + 1] = (d[i + 1] - d[i] - c[i] * h) / hi1;
        }
    }

    for (i = 0; i < n - 1; i++) {
        h = x[i + 1] - x[i];
        bi = d[i] - c[i] * h;

        a[3 * i + 0] = f[i];
        a[3 * i + 1] = bi;
        a[3 * i + 2] = c[i];
    }
}

double eval_approx2(double t, double xl, double xr, int n, const double *x, const double *a)
{
    int lo, hi, mid, i;
    double dx;

    (void)xl;
    (void)xr;

    if (n < 2) {
        return 0.0;
    }

    if (t <= x[0]) {
        i = 0;
    } else if (t >= x[n - 1]) {
        i = n - 2;
    } else {
        lo = 0;
        hi = n - 2;
        while (lo < hi) {
            mid = lo + (hi - lo) / 2;
            if (x[mid + 1] < t) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        i = lo;
    }

    dx = t - x[i];
    return a[3 * i + 0] + dx * (a[3 * i + 1] + dx * a[3 * i + 2]);
}

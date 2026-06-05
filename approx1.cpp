#include "approx1.h"
#include <math.h>


static double akima_slope(double dm2, double dm1, double d0, double d1)
{
    double w1 = fabs(d1 - d0);
    double w2 = fabs(dm2 - dm1);

    if (w1 + w2 > 1e-15 * (fabs(d0) + fabs(dm1) + 1.0)) {
        return (w1 * dm1 + w2 * d0) / (w1 + w2);
    }
    return 0.5 * (dm1 + d0);
}

void build_approx1(int n, const double *x, const double *f, double *a, double *work)
{
    int i;
    double h, di;
    double *d = work;       /* differences, size n-1 */
    double *s = work + n;   /* slopes, size n */

    if (n < 2) {
        return;
    }

    /* Step 1: compute differences */
    for (i = 0; i < n - 1; i++) {
        h = x[i + 1] - x[i];
        d[i] = (f[i + 1] - f[i]) / h;
    }

    if (n == 2) {
        s[0] = d[0];
        s[1] = d[0];
    } else if (n == 3) {
        /* Only one interior node i=1, use average */
        s[1] = 0.5 * (d[0] + d[1]);
        s[0] = (3.0 * d[0] - s[1]) / 2.0;
        s[2] = (3.0 * d[1] - s[1]) / 2.0;
    } else {
        /* Step 2: interior Akima slopes i=2..n-3 */
        for (i = 2; i <= n - 3; i++) {
            s[i] = akima_slope(d[i - 2], d[i - 1], d[i], d[i + 1]);
        }

        /* Step 3: near-boundary nodes with one ghost difference */
        {
            /* i=1: ghost d[-1] = 2*d[0] - d[1] */
            double dm1 = 2.0 * d[0] - d[1];

            s[1] = akima_slope(dm1, d[0], d[1], d[2]);
        }
        {
            /* i=n-2: ghost d[n-1] = 2*d[n-2] - d[n-3] */
            double dn1 = 2.0 * d[n - 2] - d[n - 3];

            s[n - 2] = akima_slope(d[n - 4], d[n - 3], d[n - 2], dn1);
        }

        /* Step 4: natural BC at endpoints */
        s[0] = (3.0 * d[0] - s[1]) / 2.0;
        s[n - 1] = (3.0 * d[n - 2] - s[n - 2]) / 2.0;
    }

    /* Step 5: build cubic Hermite coefficients */
    for (i = 0; i < n - 1; i++) {
        h = x[i + 1] - x[i];
        di = (f[i + 1] - f[i]) / h;

        a[4 * i + 0] = f[i];
        a[4 * i + 1] = s[i];
        a[4 * i + 2] = (3.0 * di - 2.0 * s[i] - s[i + 1]) / h;
        a[4 * i + 3] = (s[i] + s[i + 1] - 2.0 * di) / (h * h);
    }
}

double eval_approx1(double t, double xl, double xr, int n, const double *x, const double *a)
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
    return a[4 * i + 0] +
           dx * (a[4 * i + 1] + dx * (a[4 * i + 2] + dx * a[4 * i + 3]));
}

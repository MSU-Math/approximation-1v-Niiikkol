#include "approx2.h"
#include <stdlib.h>

static void solve_tridiagonal(int n, double *sub, double *diag, double *sup, double *rhs, double *x) {
    double *alpha = (double*)malloc((n-1) * sizeof(double));
    double *beta = (double*)malloc((n-1) * sizeof(double));
    
    if (!alpha || !beta) {
        free(alpha); free(beta);
        return;
    }
    
    alpha[0] = -sup[0] / diag[0];
    beta[0] = rhs[0] / diag[0];
    
    for (int i = 1; i < n-1; i++) {
        double denom = diag[i] + sub[i-1] * alpha[i-1];
        alpha[i] = -sup[i] / denom;
        beta[i] = (rhs[i] - sub[i-1] * beta[i-1]) / denom;
    }
    
    x[n-1] = (rhs[n-1] - sub[n-2] * beta[n-2]) / (diag[n-1] + sub[n-2] * alpha[n-2]);
    
    for (int i = n-2; i >= 0; i--) {
        x[i] = alpha[i] * x[i+1] + beta[i];
    }
    
    free(alpha);
    free(beta);
}

void build_approx2(int n, const double *x, const double *f, double *a, double *work) {
    (void)work;
    
    if (n < 2) return;
    
    double *d = (double*)malloc(n * sizeof(double));
    double *diag = (double*)malloc(n * sizeof(double));
    double *sub = (double*)malloc((n-1) * sizeof(double));
    double *sup = (double*)malloc((n-1) * sizeof(double));
    double *rhs = (double*)malloc(n * sizeof(double));
    double *div_diff = (double*)malloc((n-1) * sizeof(double));
    
    if (!d || !diag || !sub || !sup || !rhs || !div_diff) {
        free(d); free(diag); free(sub); free(sup); free(rhs); free(div_diff);
        return;
    }
    
    for (int i = 0; i < n-1; i++) {
        div_diff[i] = (f[i+1] - f[i]) / (x[i+1] - x[i]);
    }
    
    if (n == 2) {
        a[0] = f[0];
        a[1] = div_diff[0];
        a[2] = 0.0;
        free(d); free(diag); free(sub); free(sup); free(rhs); free(div_diff);
        return;
    }
    
    for (int i = 0; i < n; i++) {
        diag[i] = 0.0;
        rhs[i] = 0.0;
        if (i < n-1) {
            sub[i] = 0.0;
            sup[i] = 0.0;
        }
    }
    
    for (int i = 1; i < n-1; i++) {
        double h_left = x[i] - x[i-1];
        double h_right = x[i+1] - x[i];
        
        sub[i-1] = h_right;
        diag[i] = 2.0 * (h_left + h_right);
        sup[i] = h_left;
        
        rhs[i] = 3.0 * div_diff[i-1] * h_right + 3.0 * div_diff[i] * h_left;
    }
    
    diag[0] = 1.0;
    sup[0] = 1.0;
    rhs[0] = 2.0 * div_diff[0];
    
    sub[n-2] = 1.0;
    diag[n-1] = 1.0;
    rhs[n-1] = 2.0 * div_diff[n-2];
    
    solve_tridiagonal(n, sub, diag, sup, rhs, d);
    
    for (int i = 0; i < n-1; i++) {
        double h = x[i+1] - x[i];
        double a0 = f[i];
        double a1 = d[i];
        double a2 = (f[i+1] - f[i] - d[i] * h) / (h * h);
        
        a[3*i + 0] = a0;
        a[3*i + 1] = a1;
        a[3*i + 2] = a2;
    }
    
    free(d);
    free(diag);
    free(sub);
    free(sup);
    free(rhs);
    free(div_diff);
}

double eval_approx2(double t, double, double, int n, const double *x, const double *a) {
    if (n < 2) return 0.0;
    
    int i;
    if (t <= x[0]) {
        i = 0;
    } else if (t >= x[n-1]) {
        i = n-2;
    } else {
        int lo = 0, hi = n-2;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (x[mid] <= t) {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        i = lo;
    }
    
    double dt = t - x[i];
    return a[3*i + 0] + dt * (a[3*i + 1] + dt * a[3*i + 2]);
}
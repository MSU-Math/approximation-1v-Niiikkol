#include "approx2.h"
#include <stdlib.h>
#include <math.h>

static void solve_tridiagonal(int n, double *a, double *b, double *c, double *d, double *x) {
    double *alpha = (double*)malloc((n-1) * sizeof(double));
    double *beta  = (double*)malloc((n-1) * sizeof(double));
    
    alpha[0] = -c[0] / b[0];
    beta[0]  =  d[0] / b[0];
    
    for (int i = 1; i < n-1; i++) {
        double denom = b[i] + a[i] * alpha[i-1];
        alpha[i] = -c[i] / denom;
        beta[i]  = (d[i] - a[i] * beta[i-1]) / denom;
    }
    
    x[n-1] = (d[n-1] - a[n-1] * beta[n-2]) / (b[n-1] + a[n-1] * alpha[n-2]);
    
    for (int i = n-2; i >= 0; i--) {
        x[i] = alpha[i] * x[i+1] + beta[i];
    }
    
    free(alpha);
    free(beta);
}

void build_approx2(int n, const double *x, const double *f, double *a, double *work) {
    if (n < 2) return;

    double *ksi = work;        
    double *v   = work + (n+2); 
    
    ksi[0] = x[0];   
    ksi[n+1] = x[n-1]; 
    for (int i = 1; i <= n; i++) {
        ksi[i] = (x[i-1] + x[i]) / 2.0; 
    }

    ksi[0] = x[0];
    for (int i = 1; i <= n; i++) {
        ksi[i] = (x[i-1] + x[i]) / 2.0;
    }
    ksi[n] = x[n-1]; 

    int m = n + 1;
    double *diag = (double*)calloc(m, sizeof(double));
    double *sub  = (double*)calloc(m, sizeof(double));
    double *sup  = (double*)calloc(m, sizeof(double));
    double *rhs  = (double*)calloc(m, sizeof(double));
    
    for (int i = 2; i <= n; i++) {
        int idx = i - 1;  
        
        double xi_prev = x[i-2];   
        double xi_curr = x[i-1];   
        
        double ksi_prev = ksi[i-2]; 
        double ksi_curr = ksi[i-1]; 
        double ksi_next = ksi[i];   
        

        double coeff_left   = 1.0/(xi_prev - ksi_prev) - 1.0/(ksi_curr - ksi_prev);
        double coeff_center = 1.0/(ksi_curr - xi_prev) + 1.0/(ksi_curr - ksi_prev) 
                            + 1.0/(xi_curr - ksi_curr) + 1.0/(ksi_next - ksi_curr);
        double coeff_right  = 1.0/(ksi_next - xi_curr) - 1.0/(ksi_next - ksi_curr);
        
        sub[idx-1] += coeff_left;
        diag[idx]   += coeff_center;
        sup[idx+1]  += coeff_right;
        
        rhs[idx] = (1.0/(xi_prev - ksi_prev) + 1.0/(ksi_curr - xi_prev)) * f[i-2]
                 + (1.0/(xi_curr - ksi_curr) + 1.0/(ksi_next - xi_curr)) * f[i-1];
    }
    

    if (n >= 3) {
        double xi1 = x[0], xi2 = x[1];
        double ksi1 = ksi[0], ksi2 = ksi[1], ksi3 = ksi[2];
        
        double denom1 = ksi2 - ksi1;
        double denom2 = ksi3 - ksi2;
        

        double coeff_v1 = -1.0/(denom1 * (xi1 - ksi1));
        double coeff_v2 =  1.0/(denom1 * (xi1 - ksi1)) 
                         + 1.0/(denom2 * (xi2 - ksi2));
        double coeff_v3 = -1.0/(denom2 * (xi2 - ksi2));
        
        double rhs_val = -1.0/(denom1 * (ksi2 - xi1)) * f[0]
                        + (1.0/(denom1 * (ksi2 - xi1)) - 1.0/(denom2 * (ksi3 - xi2))) * f[1]
                        + 1.0/(denom2 * (ksi3 - xi2)) * f[2];
        
        sub[0] += coeff_v1;
        diag[1] += coeff_v2;
        sup[2] += coeff_v3;
        rhs[1] += rhs_val;
    }

    if (n >= 3) {
        int i = n - 2;
        double xi_prev = x[i-1], xi_curr = x[i];
        double ksi_prev = ksi[i-1], ksi_curr = ksi[i], ksi_next = ksi[i+1];
        
        double denom1 = ksi_curr - ksi_prev;
        double denom2 = ksi_next - ksi_curr;
        
        double coeff_v_prev = -1.0/(denom1 * (xi_prev - ksi_prev));
        double coeff_v_curr =  1.0/(denom1 * (xi_prev - ksi_prev)) 
                             + 1.0/(denom2 * (xi_curr - ksi_curr));
        double coeff_v_next = -1.0/(denom2 * (xi_curr - ksi_curr));
        
        double rhs_val = -1.0/(denom1 * (ksi_curr - xi_prev)) * f[i-1]
                        + (1.0/(denom1 * (ksi_curr - xi_prev)) - 1.0/(denom2 * (ksi_next - xi_curr))) * f[i]
                        + 1.0/(denom2 * (ksi_next - xi_curr)) * f[i+1];
        
        int idx = i;
        sub[idx-1] += coeff_v_prev;
        diag[idx]   += coeff_v_curr;
        sup[idx+1]  += coeff_v_next;
        rhs[idx]    += rhs_val;
    }
    

    solve_tridiagonal(m, sub, diag, sup, rhs, v + 1);

    v[0] = f[0];
    v[n] = f[n-1];
    
  
    for (int i = 0; i < n-1; i++) {
        double xi = x[i];
        double xip1 = x[i+1];
        double h = xip1 - xi;

        double x_mid = (xi + xip1) / 2.0;
        

        int k = 0;
        for (k = 0; k < n; k++) {
            if (x_mid >= ksi[k] && x_mid <= ksi[k+1]) break;
        }

        double vL = v[k];
        double vR = v[k+1];
        double xk = x[k];
        double ksiL = ksi[k];
        double ksiR = ksi[k+1];
        
        double A = (f[k] - vL) / (xk - ksiL);
        double B = 1.0/(ksiR - ksiL) * ((vR - f[k])/(ksiR - xk) - (f[k] - vL)/(xk - ksiL));
        double P_at_mid = vL + A*(x_mid - ksiL) + B*(x_mid - ksiL)*(x_mid - xk);

        double c0 = f[i];
        double c1 = (f[i+1] - f[i]) / h;
        double c2 = (P_at_mid - c0 - c1 * (x_mid - xi)) / ((x_mid - xi) * (x_mid - xi));
        
        a[3*i + 0] = c0;
        a[3*i + 1] = c1;
        a[3*i + 2] = c2;
    }
    
    free(diag);
    free(sub);
    free(sup);
    free(rhs);
}

double eval_approx2(double t, double /*xl*/, double /*xr*/, int n, const double *x, const double *a) {
    if (n < 2) return 0.0;
    
    int i;
    if (t <= x[0]) {
        i = 0;
    } else if (t >= x[n-1]) {
        i = n-2;
    } else {
        int lo = 0, hi = n-2;
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            if (x[mid+1] < t) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        i = lo;
    }
    
    double dt = t - x[i];
    return a[3*i + 0] + dt * (a[3*i + 1] + dt * a[3*i + 2]);
}
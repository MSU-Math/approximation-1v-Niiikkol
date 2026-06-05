#include "approx2.h"
#include <stdlib.h>

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
    

    int m = n+1;
    double *diag = (double*)malloc(m * sizeof(double));
    double *sub  = (double*)malloc(m * sizeof(double));
    double *sup  = (double*)malloc(m * sizeof(double));
    double *rhs  = (double*)malloc(m * sizeof(double));
    
    for (int i = 0; i < m; i++) {
        diag[i] = 0.0;
        sub[i] = 0.0;
        sup[i] = 0.0;
        rhs[i] = 0.0;
    }
    
    
    for (int i = 1; i <= n-1; i++) {
        double hL = ksi[i] - ksi[i-1];
        double hR = ksi[i+1] - ksi[i];
        double xL = x[i-1];
        double xR = x[i];
        
        double coeff_v_left   = 1.0/(xL - ksi[i-1]) - 1.0/hL;
        double coeff_v_center = 1.0/(ksi[i] - xL) + 1.0/hL + 1.0/(xR - ksi[i]) + 1.0/hR;
        double coeff_v_right  = 1.0/(ksi[i+1] - xR) - 1.0/hR;
        
        sub[i-1] = coeff_v_left;
        diag[i]  = coeff_v_center;
        sup[i+1] = coeff_v_right;
        
        rhs[i] = (1.0/(xL - ksi[i-1]) + 1.0/(ksi[i] - xL)) * f[i-1]
               + (1.0/(xR - ksi[i]) + 1.0/(ksi[i+1] - xR)) * f[i];
    }
    
    if (n >= 3) {
        double h12 = ksi[2] - ksi[1];
        double h01 = ksi[1] - ksi[0];
        double x0 = x[0], x1 = x[1];
        
        double left_v0 = -1.0/(h01 * (x0 - ksi[0]));
        double left_v1 =  1.0/(h01 * (x0 - ksi[0])) + 1.0/(h12 * (x1 - ksi[1])) - 1.0/(h12 * (ksi[2] - x1));
        double left_v2 = -1.0/(h12 * (x1 - ksi[1]));
        
        sub[0] += left_v0;
        diag[1] += left_v1;
        sup[2] += left_v2;

    }

    if (n >= 3) {
        int n_idx = n-1;
        double h_prev = ksi[n_idx] - ksi[n_idx-1];
        double h_last = ksi[n_idx+1] - ksi[n_idx];
        double x_last1 = x[n-2], x_last2 = x[n-1];
        
        double right_v_nm1 = -1.0/(h_prev * (x_last1 - ksi[n_idx-1]));
        double right_v_n   =  1.0/(h_prev * (x_last1 - ksi[n_idx-1])) 
                           + 1.0/(h_last * (x_last2 - ksi[n_idx])) 
                           - 1.0/(h_last * (ksi[n_idx+1] - x_last2));
        double right_v_np1 = -1.0/(h_last * (x_last2 - ksi[n_idx]));
        
        sub[n_idx-1] += right_v_nm1;
        diag[n_idx]   += right_v_n;
        sup[n_idx+1]  += right_v_np1;
    }
    

    solve_tridiagonal(m, sub, diag, sup, rhs, v);

    for (int i = 0; i < n-1; i++) {
        double xi = x[i];
        double xip1 = x[i+1];
        double h = xip1 - xi;
        double mid = (xi + xip1) / 2.0;

        double P_i_mid, P_ip1_mid;
        
        double ksiL = ksi[i];
        double ksiR = ksi[i+1];
        double x_center = x[i];
        double vL = v[i];
        double vR = v[i+1];
        double A = (f[i] - vL) / (x_center - ksiL);
        double B = (1.0/(ksiR - ksiL)) * ((vR - f[i])/(ksiR - x_center) - (f[i] - vL)/(x_center - ksiL));
        P_i_mid = vL + A*(mid - ksiL) + B*(mid - ksiL)*(mid - x_center);
        
        if (i < n-2) {
            ksiL = ksi[i+1];
            ksiR = ksi[i+2];
            x_center = x[i+1];
            vL = v[i+1];
            vR = v[i+2];
            A = (f[i+1] - vL) / (x_center - ksiL);
            B = (1.0/(ksiR - ksiL)) * ((vR - f[i+1])/(ksiR - x_center) - (f[i+1] - vL)/(x_center - ksiL));
            P_ip1_mid = vL + A*(mid - ksiL) + B*(mid - ksiL)*(mid - x_center);
        } else {
            P_ip1_mid = f[i+1];
        }
        
        double f_mid = (P_i_mid + P_ip1_mid) / 2.0;
        

        double c0 = f[i];
        double c1 = (f[i+1] - f[i]) / h;
        double c2 = (f_mid - c0 - c1 * h / 2.0) / (h * h / 4.0);
        
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
    int i;
    if (n < 2) {
        return 0.0;
    }
    
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
#include "poissonSolver.h"
#include <stdio.h>
 
extern double phi[GRID_SIZE][GRID_SIZE][GRID_SIZE], delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];
 
inline static double greenFunc(int,int,int,double);
 
void poissonSolver(double a) {
    const unsigned int fft_size = GRID_SIZE*GRID_SIZE*(GRID_SIZE/2+1);
    const unsigned int real_size = GRID_SIZE*GRID_SIZE*GRID_SIZE;
    fftw_complex *F_delta = (fftw_complex*)fftw_malloc(fft_size*sizeof(fftw_complex));
    fftw_plan r2c = fftw_plan_dft_r2c_3d(GRID_SIZE,GRID_SIZE,GRID_SIZE,
                    &delta[0][0][0], F_delta, FFTW_ESTIMATE);
    fftw_plan c2r = fftw_plan_dft_c2r_3d(GRID_SIZE,GRID_SIZE,GRID_SIZE,
            F_delta, &phi[0][0][0], FFTW_ESTIMATE);
 
    fftw_execute(r2c);
 
    int i,j,k;
    for (i=0; i<GRID_SIZE; i++) {
        for (j=0; j<GRID_SIZE; j++) {
            for (k=0; k<GRID_SIZE/2+1; k++) {
                register double tmp = greenFunc(i,j,k,a);
                register int index = (i*GRID_SIZE + j)*(GRID_SIZE/2+1) + k;
                F_delta[index][0] *= tmp;
                F_delta[index][1] *= tmp;
            }
        }
    }
 
    fftw_execute(c2r);
 
    double *p;
    for (i=0, p=&phi[0][0][0]; i<real_size; i++,p++) {
        *p /= real_size;
    }
 
    fftw_free(F_delta);
    fftw_destroy_plan(r2c);
    fftw_destroy_plan(c2r);
}
 
inline static double greenFunc(int l, int m, int n, double a) {
    double s_x = sin(M_PI*l/GRID_SIZE);
    double s_y = sin(M_PI*m/GRID_SIZE);
    double s_z = sin(M_PI*n/GRID_SIZE);
 
    if (l == 0 && m == 0 && n == 0) return 0.0f;
    return -3./(8.*a)/(s_x*s_x + s_y*s_y + s_z*s_z);
}
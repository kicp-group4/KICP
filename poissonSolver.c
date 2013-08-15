#include "poissonSolver.h"
 
extern double phi[GRID_SIZE][GRID_SIZE][GRID_SIZE], delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

enum {
	FFT_SIZE = GRID_SIZE*GRID_SIZE*(GRID_SIZE/2+1),
	REAL_SIZE = GRID_SIZE*GRID_SIZE*GRID_SIZE
};
static fftw_complex *F_delta;
static fftw_plan r2c, c2r;
static double* greenFunc;
 
void poissonSolver(double a) {
    fftw_execute(r2c);

    int i,j,k;
//	#pragma omp parallel for private(i,j,k)
    for(i=0;i<GRID_SIZE;i++){
    	for(j=0;j<GRID_SIZE;j++){
    		for(k=0;k<(GRID_SIZE/2+1);k++){
    			register const int index = (i*GRID_SIZE + j)*(GRID_SIZE/2+1) + k;
    			F_delta[index][0] *= greenFunc[index] / a;
    			F_delta[index][1] *= greenFunc[index] / a;
			}
		}
    }

//    int i;
//	#pragma omp parallel for private(i)
//    for(i=0; i<FFT_SIZE; i++){
//    	register double tmp = greenFunc[i]/a;
//		F_delta[i][0] *= tmp;
//		F_delta[i][1] *= tmp;
//    }
 
    fftw_execute(c2r);
 
    double *p = &phi[0][0][0];
	#pragma omp parallel for private (i)
    for (i=0; i<REAL_SIZE; i+=4) {
        p[i] /= REAL_SIZE; p[i+1] /= REAL_SIZE;
        p[i+2] /= REAL_SIZE; p[i+3] /= REAL_SIZE;
    }
}

void poissonSolver_init() {
	F_delta = (fftw_complex*)fftw_malloc(FFT_SIZE*sizeof(fftw_complex));
	greenFunc = (double*)fftw_malloc(FFT_SIZE*sizeof(double));
	r2c = fftw_plan_dft_r2c_3d(GRID_SIZE,GRID_SIZE,GRID_SIZE, &delta[0][0][0], F_delta, FFTW_ESTIMATE);
    c2r = fftw_plan_dft_c2r_3d(GRID_SIZE,GRID_SIZE,GRID_SIZE, F_delta, &phi[0][0][0], FFTW_ESTIMATE);

	int l,m,n;
	for (l=0; l<GRID_SIZE; l++) {
		for (m=0; m<GRID_SIZE; m++) {
			for (n=0; n<GRID_SIZE/2+1; n++) {
				register int index = (l*GRID_SIZE + m)*(GRID_SIZE/2+1) + n;
				if(index == 0){
					greenFunc[0] = 0.0;
					continue;
				}
				double s_x = sin(M_PI*l/GRID_SIZE);
				double s_y = sin(M_PI*m/GRID_SIZE);
				double s_z = sin(M_PI*n/GRID_SIZE);
				greenFunc[index] = -3*cosmology->OmegaM/(8.)/(s_x*s_x + s_y*s_y + s_z*s_z);
			}
		}
	}
}

void poissonSolver_cleanup() {
	fftw_free(F_delta);
	fftw_free(greenFunc);
	fftw_destroy_plan(r2c);
	fftw_destroy_plan(c2r);
}

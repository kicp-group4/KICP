#include "ic.h"
#include "cosmo_zeldovich.h"

extern struct vec3D pos, momentum;

double lambda_ijk[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double disp_ijk_x[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double disp_ijk_y[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double disp_ijk_z[GRID_SIZE][GRID_SIZE][GRID_SIZE];
static double *Pk, *k;

static const unsigned int fft_size = GRID_SIZE * GRID_SIZE * (GRID_SIZE / 2 + 1);
static const unsigned int real_size = GRID_SIZE * GRID_SIZE * GRID_SIZE;
static fftw_complex *F_lambda;
static fftw_complex *F_disp_x;
static fftw_complex *F_disp_y;
static fftw_complex *F_disp_z;
static fftw_complex *F_delta;
static fftw_plan re2co, co2re;

void cosmo_zeldovich(double a) {
	int i, j, m, n, l, n_k;
	double D_dot, D_plus, Om_Mz, Om_Lz, z;
	char file_line[50];
	FILE *input;

	z = 1.0 / a - 1.0;
	
	Om_Mz = OMEGA_M * (1.+z)*(1.+z)*(1.+z) / (OMEGA_L + OMEGA_M*(1.+z)+(1.+z)*(1.+z));
	Om_Lz = OMEGA_L / (OMEGA_L + OMEGA_M*(1.+z)+(1.+z)*(1.+z));
	
	D_plus = 5.0 / 2.0 * Om_Mz / ( pow(Om_Mz,4./7.) - Om_Lz + (1. + Om_Mz / 2) * (1. + Om_Lz / 70.) );
	D_dot = (H_0 * sqrt(Om_Mz + Om_Lz * pow(a, 3.0)) / sqrt(a));

	printf("%g\n",D_plus);
	printf("%g\n",D_dot);
	
	gsl_rng *rng = gsl_rng_alloc( gsl_rng_mt19937 );
	for (i = 0; i < N_P_1D; i++) {
		for (j = 0; j < N_P_1D; j++) {
			for (m = 0; m < N_P_1D; m++) {
				lambda_ijk[i][j][m] = gsl_ran_gaussian(rng,1.0);
			}
		}
	}
	gsl_rng_free(rng);
	
	// READING INPUT P(K)
	n_k = 0;
	input = fopen("PK_IC", "r");
	while(fgets(file_line,50,input) != NULL){
	  n_k++;
	}
	rewind(input);
	
	cosmo_zeldovich_init(n_k);

	for(j=0;j<n_k;j++){
		fscanf(input, "%lf\t%lf\n", &k[j], &Pk[j]);
		k[i] /= 1.;
		Pk[i] /= 1.;
	}
	fclose(input);

	// FFT LAMBDA
	re2co = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, &lambda_ijk[0][0][0], F_lambda, FFTW_ESTIMATE);
	fftw_execute(re2co);

	register int index;
	register const double delta_k = 2 * M_PI / GRID_SIZE;
	register double factor, tmp;
	register double imag, real;
	for (l = 0; l < GRID_SIZE; l++) {
		for (m = 0; m < GRID_SIZE; m++) {
			for (n = 0; n < GRID_SIZE / 2 + 1; n++) {
				index = (l * GRID_SIZE + m) * (GRID_SIZE / 2 + 1) + n;
				i = 0;
				tmp = delta_k * sqrt(l * l + m * m + n * n);
				while(tmp > k[i]) i++;
				if(i < n_k && i >= 0){
					F_delta[index][0] = sqrt(Pk[i]) * F_lambda[index][0]; 
					F_delta[index][1] = sqrt(Pk[i]) * F_lambda[index][1];
					
					factor = delta_k / D_plus / (tmp*tmp);
					if(l == 0 && m == 0 && n == 0) factor = 0.;
					real = factor * F_delta[index][0];
					imag = factor * F_delta[index][1];
					
					F_disp_x[index][1] = - l * real;
					F_disp_y[index][1] = - m * real;
					F_disp_z[index][1] = - n * real;

					F_disp_x[index][0] = l * imag;
					F_disp_y[index][0] = m * imag;
					F_disp_z[index][0] = n * imag;
				}
			}
		}
	}

	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_x, &disp_ijk_x[0][0][0], FFTW_ESTIMATE);
	fftw_execute(co2re);
	
	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_y, &disp_ijk_y[0][0][0], FFTW_ESTIMATE);
	fftw_execute(co2re);
	
	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_z, &disp_ijk_z[0][0][0], FFTW_ESTIMATE);
	fftw_execute(co2re);
	
	double *px,*py,*pz;
	for (i=0, px=&disp_ijk_x[0][0][0], py=&disp_ijk_y[0][0][0], pz=&disp_ijk_z[0][0][0]; i<real_size; i++,px++,py++,pz++) {
	    *px /= real_size;
    	    *py /= real_size;
	    *pz /= real_size;
	}

	for (i = 0; i < N_P_1D; i++) {
		for (j = 0; j < N_P_1D; j++) {
			for (m = 0; m < N_P_1D; m++) {
				pos.x[i + N_P_1D * (j + N_P_1D * m)] = 0.5 + i + D_plus * disp_ijk_x[i][j][m];
				pos.y[i + N_P_1D * (j + N_P_1D * m)] = 0.5 + j + D_plus * disp_ijk_y[i][j][m];
				pos.z[i + N_P_1D * (j + N_P_1D * m)] = 0.5 + m + D_plus * disp_ijk_z[i][j][m];
				
				momentum.x[i + N_P_1D * (j + N_P_1D * m)] = a * D_dot * disp_ijk_x[i][j][m];
				momentum.y[i + N_P_1D * (j + N_P_1D * m)] = a * D_dot * disp_ijk_y[i][j][m];
				momentum.z[i + N_P_1D * (j + N_P_1D * m)] = a * D_dot * disp_ijk_z[i][j][m];
			}
		}
	}
	cosmo_zeldovich_cleanup();

}

void cosmo_zeldovich_init(int n_ks) {
	Pk = (double *) malloc(n_ks * sizeof(double));
	k  = (double *) malloc(n_ks * sizeof(double));
	F_lambda = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_delta  = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_x = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_y = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_z = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
}

void cosmo_zeldovich_cleanup() {
	free(k);
	free(Pk);
	fftw_free(F_delta);
	fftw_free(F_lambda);
	fftw_free(F_disp_x);
	fftw_free(F_disp_y);
	fftw_free(F_disp_z);
	fftw_destroy_plan(re2co);
	fftw_destroy_plan(co2re);
}

#include "cosmo_zeldovich.h"

extern struct vec3D pos, momentum;

static inline void cosmo_zeldovich_cleanup();
static inline void cosmo_zeldovich_init(int);

static double *lambda_ijk, *disp_ijk_x, *disp_ijk_y, *disp_ijk_z, *Pk, *k;
static const unsigned int fft_size = GRID_SIZE * GRID_SIZE * (GRID_SIZE / 2 + 1);
static const unsigned int real_size = GRID_SIZE * GRID_SIZE * GRID_SIZE;
static fftw_complex *F_lambda, *F_disp_x, *F_disp_y, *F_disp_z, *F_delta_lmn;
static fftw_plan re2co, co2re;

void cosmo_zeldovich(double a) {
	int i, j, m, n, l, n_k;
	double x, y;
	char file_line[50];
	const double scale = L_BOX / (double) GRID_SIZE;
	FILE *input;

	const double D_plus = dPlus(a);
	const double D_dot = qPlus(a) * cosmology->h * 100 / (a * a);
	const double a_temp = A_INITIAL - DELTA_A / 2.0;

	/* Can't all of this messy file business be replaced with a compiled-time constant? */
	n_k = 0;
	input = fopen("PK_IC", "r");
	while (fgets(file_line, 50, input) != NULL) {
		n_k++;
	}
	rewind(input);

	cosmo_zeldovich_init(n_k);

	/* With a compile-time constant, we can speed this up */
	for (j = 0; j < n_k; j++) {
		fscanf(input, "%lf\t%lf\n", &k[j], &Pk[j]);
		k[j] *= scale;
		Pk[j] *= scale * scale * scale;
	}
	fclose(input);

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	#pragma omp parallel for private(i)
	for(i=0; i<real_size; i+=2){
		lambda_ijk[i] = gsl_ran_gaussian(rng, 1.0);
		lambda_ijk[i+1] = gsl_ran_gaussian(rng, 1.0);
	}
	gsl_rng_free(rng);

	re2co = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, lambda_ijk, F_lambda, FFTW_ESTIMATE);
	fftw_execute(re2co);

	const double delta_k = 2 * M_PI / GRID_SIZE;

	for (l = 0, x = 0; l < GRID_SIZE; l++, x++) {
		if (l >= GRID_SIZE / 2)
			x = l - GRID_SIZE;
		for (m = 0, y = 0; m < GRID_SIZE; m++, y++) {
			if (m >= GRID_SIZE / 2)
				y = m - GRID_SIZE;
			for (n = 0; n < GRID_SIZE / 2 + 1; n++) {
				register const int index = (l * GRID_SIZE + m) * (GRID_SIZE / 2 + 1) + n;
				register const double tmp = delta_k * sqrt(x * x + y * y + n * n);
				i = 0;
				while (tmp > k[i])
					i++;
				if (i < n_k && i >= 0) {
					F_delta_lmn[index][0] = sqrt(Pk[i]) * F_lambda[index][0];
					F_delta_lmn[index][1] = sqrt(Pk[i]) * F_lambda[index][1];

					register double factor = delta_k / D_plus / (tmp * tmp);

					if (x == 0 && y == 0 && n == 0)
						factor = 0.;
					register const double real = factor * F_delta_lmn[index][0];
					register const double imag = factor * F_delta_lmn[index][1];

					F_disp_x[index][1] = -x * real;
					F_disp_y[index][1] = -y * real;
					F_disp_z[index][1] = -n * real;

					F_disp_x[index][0] = x * imag;
					F_disp_y[index][0] = y * imag;
					F_disp_z[index][0] = n * imag;
				}
			}
		}
	}

	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_x, disp_ijk_x, FFTW_ESTIMATE);
	fftw_execute(co2re);

	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_y, disp_ijk_y, FFTW_ESTIMATE);
	fftw_execute(co2re);

	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_z, disp_ijk_z, FFTW_ESTIMATE);
	fftw_execute(co2re);

	#pragma omp parallel for private(i)
	for(i=0; i<real_size; i+=2){
		disp_ijk_x[i] /= real_size; disp_ijk_x[i+1] /= real_size;
		disp_ijk_y[i] /= real_size; disp_ijk_y[i+1] /= real_size;
		disp_ijk_z[i] /= real_size; disp_ijk_z[i+1] /= real_size;
	}

	/* --- DANGER --- : What if N_P_1D != GRID_SIZE ?  */

	register double bob = GRID_SIZE / (double) N_P_1D;
	for (i = 0; i < N_P_1D; i++) {
		for (j = 0; j < N_P_1D; j++) {
			for (m = 0; m < N_P_1D; m++) {
				register const int p_index = i + N_P_1D * (j + N_P_1D * m);			/* for pos/momentum */
				register const int d_index = i + GRID_SIZE * (j + GRID_SIZE * m);	/* for disp_ijk_* */
				pos.x[p_index] = bob * i + D_plus * disp_ijk_x[d_index];
				pos.y[p_index] = bob * j + D_plus * disp_ijk_y[d_index];
				pos.z[p_index] = bob * m + D_plus * disp_ijk_z[d_index];

				momentum.x[p_index] = a_temp * a_temp * D_dot * disp_ijk_x[d_index];
				momentum.y[p_index] = a_temp * a_temp * D_dot * disp_ijk_y[d_index];
				momentum.z[p_index] = a_temp * a_temp * D_dot * disp_ijk_z[d_index];

				while (pos.x[p_index] > (double) GRID_SIZE)
					pos.x[p_index] -= (double) GRID_SIZE;
				while (pos.y[p_index] > (double) GRID_SIZE)
					pos.y[p_index] -= (double) GRID_SIZE;
				while (pos.z[p_index] > (double) GRID_SIZE)
					pos.z[p_index] -= (double) GRID_SIZE;
				while (pos.x[p_index] < 0.)
					pos.x[p_index] += (double) GRID_SIZE;
				while (pos.y[p_index] < 0.)
					pos.y[p_index] += (double) GRID_SIZE;
				while (pos.z[p_index] < 0.)
					pos.z[p_index] += (double) GRID_SIZE;
			}
		}
	}

	cosmo_zeldovich_cleanup();
}

static inline void cosmo_zeldovich_init(int n_ks) {
	F_lambda = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_delta_lmn = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_x = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_y = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_z = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));

	/* Use fftw_malloc for the 'real' arrays to get good alignment */
	Pk = (double *) fftw_malloc(n_ks * sizeof(double));
	k = (double *) fftw_malloc(n_ks * sizeof(double));
	lambda_ijk = (double*) fftw_malloc(real_size * sizeof(double));
	disp_ijk_x = (double*) fftw_malloc(real_size * sizeof(double));
	disp_ijk_y = (double*) fftw_malloc(real_size * sizeof(double));
	disp_ijk_z = (double*) fftw_malloc(real_size * sizeof(double));
}

static inline void cosmo_zeldovich_cleanup() {
	fftw_free(k);
	fftw_free(Pk);
	fftw_free(F_delta_lmn);
	fftw_free(F_lambda);
	fftw_free(F_disp_x);
	fftw_free(F_disp_y);
	fftw_free(F_disp_z);
	fftw_free(lambda_ijk);
	fftw_free(disp_ijk_x);
	fftw_free(disp_ijk_y);
	fftw_free(disp_ijk_z);
	fftw_destroy_plan(re2co);
	fftw_destroy_plan(co2re);
}

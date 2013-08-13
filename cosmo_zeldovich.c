#include "ic.h"
#include "cosmo_zeldovich.h"

extern struct vec3D pos, momentum;

static inline void cosmo_zeldovich_cleanup();
static inline void cosmo_zeldovich_init(int);

double lambda_ijk[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double disp_ijk_x[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double disp_ijk_y[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double disp_ijk_z[GRID_SIZE][GRID_SIZE][GRID_SIZE];
static double *Pk, *k;

static const unsigned int fft_size = GRID_SIZE * GRID_SIZE
		* (GRID_SIZE / 2 + 1);
static const unsigned int real_size = GRID_SIZE * GRID_SIZE * GRID_SIZE;
static fftw_complex *F_lambda;
static fftw_complex *F_disp_x;
static fftw_complex *F_disp_y;
static fftw_complex *F_disp_z;
static fftw_complex *F_delta_lmn;
static fftw_plan re2co, co2re;

void cosmo_zeldovich(double a) {
	int i, j, m, n, l, n_k;
	double x, y;
	double D_dot, D_plus;
	char file_line[50];
	double scale = L_BOX / (double) GRID_SIZE;
	FILE *input;

	cosmology_set(OmegaM, OMEGA_M);
	cosmology_set(OmegaL, OMEGA_L);
	cosmology_set(OmegaB, 0.04);
	cosmology_set(h, 0.702);
	cosmology_set_thread_safe_range(1.0e-3, 1);

	D_plus = dPlus(a);
	double a_temp = A_INITIAL - DELTA_A / 2.0;

	D_dot = sqrt(OMEGA_M + OMEGA_L * pow(a_temp, 3.0)) / sqrt(a_temp);

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	for (i = 0; i < GRID_SIZE; i++) {
		for (j = 0; j < GRID_SIZE; j++) {
			for (m = 0; m < GRID_SIZE; m++) {
				lambda_ijk[i][j][m] = gsl_ran_gaussian(rng, 1.0);
			}
		}
	}
	gsl_rng_free(rng);

	n_k = 0;
	input = fopen("PK_IC", "r");
	while (fgets(file_line, 50, input) != NULL) {
		n_k++;
	}
	rewind(input);

	cosmo_zeldovich_init(n_k);
	for (j = 0; j < n_k; j++) {
		fscanf(input, "%lf\t%lf\n", &k[j], &Pk[j]);
		k[j] *= scale;
		Pk[j] *= scale * scale * scale;
	}
	fclose(input);

	re2co = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE,
			&lambda_ijk[0][0][0], F_lambda, FFTW_ESTIMATE);
	fftw_execute(re2co);

	register int index;
	register const double delta_k = 2 * M_PI / GRID_SIZE;
	register double factor, tmp;
	register double imag, real;

	for (l = 0, x = 0; l < GRID_SIZE; l++, x++) {

		if (l >= GRID_SIZE / 2)
			x = l - GRID_SIZE;
		for (m = 0, y = 0; m < GRID_SIZE; m++, y++) {

			if (m >= GRID_SIZE / 2)
				y = m - GRID_SIZE;
			for (n = 0; n < GRID_SIZE / 2 + 1; n++) {
				index = (l * GRID_SIZE + m) * (GRID_SIZE / 2 + 1) + n;

				tmp = delta_k * sqrt(x * x + y * y + n * n);
				i = 0;
				while (tmp > k[i])
					i++;
				if (i < n_k && i >= 0) {
					F_delta_lmn[index][0] = sqrt(Pk[i]) * F_lambda[index][0];
					F_delta_lmn[index][1] = sqrt(Pk[i]) * F_lambda[index][1];

					factor = delta_k / D_plus / (tmp * tmp);

					if (x == 0 && y == 0 && n == 0)
						factor = 0.;
					real = factor * F_delta_lmn[index][0];
					imag = factor * F_delta_lmn[index][1];

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

	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_x,
			&disp_ijk_x[0][0][0], FFTW_ESTIMATE);
	fftw_execute(co2re);

	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_y,
			&disp_ijk_y[0][0][0], FFTW_ESTIMATE);
	fftw_execute(co2re);

	co2re = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, F_disp_z,
			&disp_ijk_z[0][0][0], FFTW_ESTIMATE);
	fftw_execute(co2re);

	double *px, *py, *pz;
	for (i = 0, px = &disp_ijk_x[0][0][0], py = &disp_ijk_y[0][0][0], pz =
			&disp_ijk_z[0][0][0]; i < real_size; i++, px++, py++, pz++) {
		*px /= real_size;
		*py /= real_size;
		*pz /= real_size;
	}
	register double bob = GRID_SIZE / (double) N_P_1D;
	for (i = 0; i < N_P_1D; i++) {
		for (j = 0; j < N_P_1D; j++) {
			for (m = 0; m < N_P_1D; m++) {
				index = i + N_P_1D * (j + N_P_1D * m);
				pos.x[index] = bob * i + D_plus * disp_ijk_x[i][j][m];
				pos.y[index] = bob * j + D_plus * disp_ijk_y[i][j][m];
				pos.z[index] = bob * m + D_plus * disp_ijk_z[i][j][m];

				momentum.x[index] = a_temp * a_temp * D_dot
						* disp_ijk_x[i][j][m];
				momentum.y[index] = a_temp * a_temp * D_dot
						* disp_ijk_y[i][j][m];
				momentum.z[index] = a_temp * a_temp * D_dot
						* disp_ijk_z[i][j][m];

				while (pos.x[index] > (double) GRID_SIZE)
					pos.x[index] -= (double) GRID_SIZE;
				while (pos.y[index] > (double) GRID_SIZE)
					pos.y[index] -= (double) GRID_SIZE;
				while (pos.z[index] > (double) GRID_SIZE)
					pos.z[index] -= (double) GRID_SIZE;
				while (pos.x[index] < 0.)
					pos.x[index] += (double) GRID_SIZE;
				while (pos.y[index] < 0.)
					pos.y[index] += (double) GRID_SIZE;
				while (pos.z[index] < 0.)
					pos.z[index] += (double) GRID_SIZE;
			}
		}
	}

	cosmo_zeldovich_cleanup();
	printf("IC DONE!\n");

}

static inline void cosmo_zeldovich_init(int n_ks) {
	Pk = (double *) malloc(n_ks * sizeof(double));
	k = (double *) malloc(n_ks * sizeof(double));
	F_lambda = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_delta_lmn = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_x = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_y = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
	F_disp_z = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
}

static inline void cosmo_zeldovich_cleanup() {
	free(k);
	free(Pk);
	fftw_free(F_delta_lmn);
	fftw_free(F_lambda);
	fftw_free(F_disp_x);
	fftw_free(F_disp_y);
	fftw_free(F_disp_z);
	fftw_destroy_plan(re2co);
	fftw_destroy_plan(co2re);
}

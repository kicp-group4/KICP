#include "power_spectrum.h"

extern double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];

static const unsigned int fft_size = GRID_SIZE * GRID_SIZE
		* (GRID_SIZE / 2 + 1);
static fftw_complex *F_rho;
static fftw_plan re2co;
static double *k_M, *Pk_M;
static int *av;

void power_spectrum() {
	FILE *output;
	int i, l, m, n;
	const double scale = L_BOX / (double) GRID_SIZE;
	const double volume = (GRID_SIZE * GRID_SIZE * GRID_SIZE);

	for (i = 0; i < NUM_BINS; i+=2) {
		Pk_M[i] = 0.0; Pk_M[i + 1] = 0.0;
		av[i] = 0; av[i + 1] = 0;
		const register double tmp1 = (M_PI - 2 * M_PI / GRID_SIZE) / NUM_BINS;
		const register double tmp2 = 2 * M_PI / GRID_SIZE;
		k_M[i] = i * tmp1 + tmp2;
		k_M[i + 1] = (i+1) * tmp1 + tmp2;
	}

	re2co = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, &rho[0][0][0], F_rho, FFTW_ESTIMATE);
	fftw_execute(re2co);

	register const double factor = 2 * M_PI / (double) GRID_SIZE;
	#pragma omp parallel for private(l,m,n,i)
	for (l = 0; l < GRID_SIZE; l++) {
		for (m = 0; m < GRID_SIZE; m++) {
			for (n = 0; n < GRID_SIZE / 2 + 1; n++) {
				register const int index = (l * GRID_SIZE + m) * (GRID_SIZE / 2 + 1) + n;
				i = 0;
				register const double tmp = factor * sqrt(l * l + m * m + n * n);
				while (tmp > k_M[i])
					i++;
				if (i < NUM_BINS && i >= 0 && l + m + n != 0) {
					Pk_M[i] += F_rho[index][0] * F_rho[index][0]
							+ F_rho[index][1] * F_rho[index][1];
					av[i]++;
				}
			}
		}
	}

	output = fopen("Pk.dat", "a");
	for (i = 0; i < NUM_BINS; i++) {
		if (av[i] != 0) {
			Pk_M[i] = Pk_M[i] / (double) av[i] / (scale * scale * scale);
			fprintf(output, "%g\t%g\n", k_M[i] / scale, Pk_M[i] / volume);
		}
	}
	fclose(output);
}

void power_spectrum_init() {
	Pk_M = (double *) malloc(NUM_BINS * sizeof(double));
	k_M = (double *) malloc(NUM_BINS * sizeof(double));
	av = (int *) malloc(NUM_BINS * sizeof(int));
	F_rho = (fftw_complex*) fftw_malloc(fft_size * sizeof(fftw_complex));
}

void power_spectrum_cleanup() {
	free(k_M);
	free(Pk_M);
	free(av);
	fftw_free(F_rho);
	fftw_destroy_plan(re2co);
}

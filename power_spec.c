#include "power_spectrum.h"

extern double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];

static const unsigned int fft_size = GRID_SIZE * GRID_SIZE * (GRID_SIZE / 2 + 1);
static fftw_complex *F_rho;
static fftw_plan re2co;
static double *k_M, *Pk_M;
static int *av;

void power_spectrum() {
	FILE *output;
	int i, l, m, n;
	double tmp;

	for (i = 0; i < NUM_BINS; i++) {
		Pk_M[i] = 0;
		k_M[i] = i+1;
		av[i] = 0;
	}

	re2co = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, &rho[0][0][0], F_rho, FFTW_ESTIMATE);
	fftw_execute(re2co);

	register const double factor = 2 * M_PI / GRID_SIZE;
	for (l = 0; l < GRID_SIZE; l++) {
		for (m = 0; m < GRID_SIZE; m++) {
			for (n = 0; n < GRID_SIZE / 2 + 1; n++) {
				register const int index = (l * GRID_SIZE + m) * (GRID_SIZE / 2 + 1) + n;
				i = 0;
				tmp = factor * sqrt(l * l + m * m + n * n);
				while (tmp > k_M[i])
					i++;
				Pk_M[i] += F_rho[index][0] * F_rho[index][0] + F_rho[index][1] * F_rho[index][1];
				av[i]++;
			}
		}
	}

	output = fopen("Pk.dat", "a");
	for (i = 0; i < NUM_BINS; i++) {
		if (av[i] == 0)
			av[i] = 1;
		Pk_M[i] /= (double) av[i];
		fprintf(output, "%lg\t%lg\t%lg\n", k_M[i], 2*M_PI*L_BOX/k_M[i], Pk_M[i]);
	}
	fclose(output);
}

void power_spectrum_init() {
	k_M = (double *) malloc(NUM_BINS * sizeof(double));
	Pk_M = (double *) malloc(NUM_BINS * sizeof(double));
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

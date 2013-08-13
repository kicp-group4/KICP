#include "constants.h"
#include "vec3D.h"
#include "output.h"
#include "updateValues.h"
#include "poissonSolver.h"
#include "power_spectrum.h"
#include "cosmo_zeldovich.h"
#include <unistd.h>
#include <sys/time.h>
#include <fftw3.h>

#ifdef ST_GLASS
	#include "ic_glass.h"
#elif defined ST_COSMO_ZELDOVICH
	#include "cosmo_zeldovich.h"
#elif defined ST_ZELDOVICH
	#include "ic.h"
#else
	# error "Must specify a simulation type"
#endif

void init();
void cleanup();

struct vec3D pos, momentum;
double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double phi[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

double a = A_INITIAL;

int main() {
	init();
	unlink("pos.dat");
	unlink("density.dat");
	unlink("Pk.dat");

#ifdef ST_ZELDOVICH
	ic(a);
#elif defined ST_GLASS
	ic_glass();
#elif defined ST_COSMO_ZELDOVICH
	cosmo_zeldovich(a);
#endif

	update_density(a);
	power_spectrum();

	struct timeval t1, t2;
	gettimeofday(&t1, 0);

	char posname[12];
	char densname[15];
	int n;
	for (n = 0; n < 100; n++) {
		poissonSolver(a);
		update_particles(a);
		update_density(a);
		if (n % 20 == 0) {
			sprintf(posname, "pos_%d.txt", n);
			sprintf(densname, "density_%d.txt", n);
			output(a, posname, densname);
		}
		a += DELTA_A;
	}
	gettimeofday(&t2, 0);
	printf("time = %lg\n", (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6);
	power_spectrum();
	cleanup();
	printf("Done!\n");
	return 0;
}

void init() {
	/************** Call these ONLY once *************/
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
	/*************************************************/
	poissonSolver_init();
	power_spectrum_init();
}

void cleanup() {
	fftw_cleanup_threads();
	poissonSolver_cleanup();
	power_spectrum_cleanup();
}

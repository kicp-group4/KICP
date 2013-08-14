#include "constants.h"
#include "vec3D.h"
#include "output.h"
#include "updateValues.h"
#include "poissonSolver.h"
#include "power_spectrum.h"
#include "cosmo_zeldovich.h"
#include "cosmology.h"
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

	char posname[20];
	int n;
	int a_step;
	for (n = 0; n < 1000; n++) {
		poissonSolver(a);
		update_particles(a);
		update_density(a);
		if (n % 50 == 0) {
		  a_step = floor(a*10000);
		  if (a_step >= 1000){
		    sprintf(posname, "pos_%d.txt", a_step);
		    output(a, posname);
		  }
		  else if (a_step >= 100){
		    sprintf(posname, "pos_0%d.txt", a_step);
		    output(a, posname);
		  }
		  else if (a_step >= 10){
		    sprintf(posname, "pos_00%d.txt", a_step);
		    output(a, posname);
		  }
		  else{
		    sprintf(posname, "pos_000%d.txt", a_step);
		    output(a, posname);
		  }
		}
		a += DELTA_A;
	}
	gettimeofday(&t2, 0);
	printf("a = %g",a);
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

#ifdef ST_COSMO_ZELDOVICH
	cosmology_set(OmegaM, 0.275);
	cosmology_set(OmegaL, 0.725);
	cosmology_set(OmegaB, 0.04);
#else
	cosmology_set(OmegaM, 1.0);
	cosmology_set(OmegaL, 0.0);
	cosmology_set(OmegaB, 0.04);
#endif

	cosmology_set(h, 0.702);
	cosmology_set_thread_safe_range(1.0e-3, 1);

	poissonSolver_init();
	power_spectrum_init();
}

void cleanup() {
	fftw_cleanup_threads();
	poissonSolver_cleanup();
	power_spectrum_cleanup();
}

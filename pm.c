#include "constants.h"
#include "vec3D.h"
#include "ic.h"
#include "output.h"
#include "updateValues.h"
#include "poissonSolver.h"
#include "power_spectrum.h"
#include <unistd.h>

#ifdef ST_GLASS
	#include "ic_glass.h"
#endif

void init();
void cleanup();

struct vec3D pos, momentum;
double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double phi[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

double a = A_INITIAL;
   
int main(){
	init();
	unlink("pos.dat");
	unlink("density.dat");
	unlink("Pk.dat");

#ifdef ST_ZELDOVICH
	ic(a);
#elif defined ST_GLASS
	ic_glass();
#else
# error "Must specify a simulation type"
#endif

	update_density(a);
	power_spectrum();

	int n;
	for(n=0;n<200;n++) {
		poissonSolver(a);
		update_particles(a);
		update_density(a);
//		output(a,"pos.dat","density.dat");
		a+=DELTA_A;
		printf("Time Step: %i\n",n);
	}
	power_spectrum();
	cleanup();
	printf("Done!\n");
    return 0;
}

void init() {
	poissonSolver_init();
	power_spectrum_init();
}

void cleanup() {
	poissonSolver_cleanup();
	power_spectrum_cleanup();
}

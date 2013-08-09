#include "constants.h"
#include "vec3D.h"
#include "ic.h"
#include "output.h"
#include "shape.h"
#include "updateValues.h"
#include "poissonSolver.h"
#include "ic_glass.h"
#include <unistd.h>

void init();
void cleanup();

struct vec3D pos, momentum;
double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double phi[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

typedef enum {ZELDOVICH,GLASS,ZELDOVICH_COSMO} SimulationType;

double a = A_INITIAL;
Shape shape = CIC;
SimulationType sim = GLASS;
double accelerationSign = 1.0;
   
int main(){
	init();
	unlink("pos.dat");
	unlink("density.dat");

	switch(sim) {
	case ZELDOVICH: ic(a); break;
	case GLASS:
		ic_glass();
		accelerationSign = -1.0;
		break;
	default:
		cleanup(); fprintf(stderr,"Bad simulation type: %d\n",sim); exit(-1);
	}

	update_density(shape,a);
	output(a,"pos.dat","density.dat");

	int n;
	for(n=0;n<20;n++) {
		poissonSolver(a);
		update_particles(shape,a,accelerationSign);
		update_density(shape,a);
		output(a,"pos.dat","density.dat");
		a+=DELTA_A;
		printf("Time Step: %i\n",n);
	}
	cleanup();
	printf("Done!\n");
    return 0;
}

void init() {
	poissonSolver_init();
}

void cleanup() {
	poissonSolver_cleanup();
}

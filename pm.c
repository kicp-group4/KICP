#include "constants.h"
#include "vec3D.h"
#include "ic.h"
#include "output.h"
#include "shape.h"
#include "updateValues.h"
#include "poissonSolver.h"

void init();
void cleanup();

struct vec3D pos, momentum;
double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double phi[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

double a = A_INITIAL;
Shape shape = CIC;
   
int main(){
	init();
	ic(a);
	printf("init done!\n");
	update_density(shape,a);
	printf("density updated!\n");
	output(a,"pos.dat","density.dat");

	int n;
	for(n=0;n<10;n++) {
		poissonSolver(a);
		printf("poisson!\n");
		
		update_particles(shape,a);
		printf("particles updated!\n");
		
		update_density(shape,a);
		printf("density updated!\n");
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

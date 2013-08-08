#include "constants.h"
#include "vec3D.h"
#include "ic.h"
#include "output.h"
#include "shape.h"
#include "updateValues.h"
#include "poissonSolver.h"

struct vec3D pos, momentum;
double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double phi[GRID_SIZE][GRID_SIZE][GRID_SIZE];
double delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

double a = A_INITIAL;
   
int main(){
	ic(a);
	printf("init done!\n");
	update_density(CIC,a);
	printf("density updated!\n");
	output(a,"pos.dat","density.dat");

	int n;
	for(n=0;n<10;n++) {
		poissonSolver(a);
		printf("poisson!\n");
		
		update_particles(CIC,a);
		printf("particles updated!\n");
		printf("works?\n");
		
		update_density(CIC,a);
		printf("density updated!\n");
		output(a,"pos.dat","density.dat");
		a+=DELTA_A;
		
		printf("Time Step \t%i\t  Done!\n",n);
	}
    return 0;
}

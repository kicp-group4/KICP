#include "ic_glass.h"

extern struct vec3D pos, momentum;

void ic_glass() {
	gsl_rng *rng = gsl_rng_alloc( gsl_rng_mt19937 );

	int i;
	for(i=0; i<N_PARTICLES;i++) {
		pos.x[i] = gsl_rng_uniform(rng)*GRID_SIZE;
		pos.y[i] = gsl_rng_uniform(rng)*GRID_SIZE;
		pos.z[i] = gsl_rng_uniform(rng)*GRID_SIZE;
		momentum.x[i] = 0;
		momentum.y[i] = 0;
		momentum.z[i] = 0;
	}
	gsl_rng_free(rng);
}

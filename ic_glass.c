#include "ic_glass.h"

extern struct vec3D pos, momentum;

void ic_glass() {
	gsl_rng *rng = gsl_rng_alloc( gsl_rng_mt19937 );

	int i,j,k;
	for(i=0; i<N_PARTICLES;i++) {
		for(j=0;j<N_PARTICLES;j++) {
			for(k=0;k<N_PARTICLES;k++) {
				pos[i][j][k] = gsl_rng_uniform(rng);
				momentum[i][j][k] = 0;
			}
		}
	}
	gsl_rng_free(rng);
}

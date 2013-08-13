#include "output.h"

extern struct vec3D pos, momentum;
extern double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];

void output(double a, char *out_xyz_file, char *out_density) {

	FILE *output;
	int i, j, k;
	static int showParameters = 1;

	// POSITION + VELOCITY 
	output = fopen(out_xyz_file, "a");

	if(showParameters) {
		fprintf(output, "%d\n%f %f %f %f %f %f \n",
				N_PARTICLES, 0.0, 0.0, 0.0,
				(float)GRID_SIZE,(float)GRID_SIZE,(float)GRID_SIZE);
	}

	for (i = 0; i < N_PARTICLES; i++) {
		fprintf(output, "%lg %lg %lg %lg %lg %lg\n",
				(float)pos.x[i] , (float)pos.y[i], (float)pos.z[i],
				(float)momentum.x[i], (float)momentum.y[i], (float)momentum.z[i]);
	}
	fclose(output);

	// DENSITY OUTPUT
	output = fopen(out_density, "a");

	if(showParameters)
		fprintf(output, "%d\n", GRID_SIZE);
	for (k = 0; k < GRID_SIZE; k++) {
		for (j = 0; j < GRID_SIZE; j++) {
			for (i = 0; i < GRID_SIZE; i++) {
				fprintf(output, "%g ", rho[i][j][k]);
			}
			fprintf(output, "\n");
		}
	}
	fclose(output);
	showParameters = 0;
}

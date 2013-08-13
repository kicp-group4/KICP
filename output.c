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
		fprintf(output, "%lg\n%d\n%d\n%lg\n%lg\n%lg\n%lg\n",
				L_BOX, GRID_SIZE, N_PARTICLES, H_0,
				OMEGA_M, OMEGA_L, M_PARTICLE);
	}

	for (i = 0; i < N_PARTICLES; i++) {
		fprintf(output, "%lg %lg %lg %lg %lg %lg\n",
				pos.x[i] , pos.y[i], pos.z[i],
				momentum.x[i], momentum.y[i], momentum.z[i]);
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

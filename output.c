#include "output.h"

extern struct vec3D pos, momentum;

void output(double a, char *out_xyz_file) {

	FILE *output;
	int i;

	output = fopen(out_xyz_file, "w");
	fprintf(output, "%d\n%f %f %f %f %f %f \n",
				N_PARTICLES, 0.0, 0.0, 0.0,
				(float) GRID_SIZE, (float) GRID_SIZE, (float) GRID_SIZE);

	for (i = 0; i < N_PARTICLES; i++) {
		fprintf(output, "%lg %lg %lg %lg %lg %lg\n",
				(float) pos.x[i], (float) pos.y[i], (float) pos.z[i],
				(float) momentum.x[i], (float) momentum.y[i], (float) momentum.z[i]);
	}
	fclose(output);
}

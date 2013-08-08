#include "output.h"

extern struct vec3D pos, momentum;
extern double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE];

void output(double a, char *out_xyz_file, char *out_density) {

	FILE *output;
	int i, j, k;
	double a_cubed = a * a * a;

	// POSITION + VELOCITY 
	output = fopen(out_xyz_file, "a");
	fprintf(output, "%g\n", L_BOX);
	fprintf(output, "%d\n", GRID_SIZE);
	fprintf(output, "%d\n", N_PARTICLES);
	fprintf(output, "%i\n", H_0);
	fprintf(output, "%g\n", OMEGA_M);
	fprintf(output, "%g\n", OMEGA_L);
	fprintf(output, "%g\n", M_PARTICLE);

	for (i = 0; i < N_PARTICLES; i++) {
		fprintf(output, "%g \t %g \t %g",
				pos.x[i] * R_0 * a/ 3.086e16,
				pos.y[i] * R_0 * a/ 3.086e16,
				pos.z[i] * R_0 * a/ 3.086e16);
		fprintf(output, " \t %g \t %g \t %g \n",
				momentum.x[i] * V_0 / a,
				momentum.y[i] * V_0 / a,
				momentum.z[i] * V_0 / a);
	}
	fprintf(output, "\n**********\n");
	fclose(output);

	// DENSITY OUTPUT
	output = fopen(out_density, "a");
	fprintf(output, "%d\n", GRID_SIZE);
	for (i = 0; i < GRID_SIZE; i++) {
		for (j = 0; j < GRID_SIZE; j++) {
			for (k = 0; k < GRID_SIZE; k++) {
				fprintf(output, "%g \t", rho[i][j][k] * RHO_0 / a_cubed);
			}
			fprintf(output, "\n");
		}
	}
	fprintf(output, "\n**********\n");
	fclose(output);

}

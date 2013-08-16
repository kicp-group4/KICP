#include "output.h"

extern struct vec3D pos;

void output(char *out_xyz_file) {
	int ntemp, n = N_PARTICLES;
	FILE *F;
	float cornerPos = 0.0;
	F = fopen(out_xyz_file, "wb");

	ntemp = sizeof(float);
	fwrite(&ntemp, 4, 1, F);
	fwrite(&n, 4, 1, F);
	fwrite(&ntemp, 4, 1, F);

	ntemp = 24;
	fwrite(&ntemp, 4, 1, F);
	fwrite(&cornerPos, 4, 1, F);
	fwrite(&cornerPos, 4, 1, F);
	fwrite(&cornerPos, 4, 1, F);
	cornerPos = (float) GRID_SIZE;
	fwrite(&cornerPos, 4, 1, F);
	fwrite(&cornerPos, 4, 1, F);
	fwrite(&cornerPos, 4, 1, F);
	fwrite(&ntemp, 4, 1, F);

	ntemp = sizeof(double) * n;
	fwrite(&ntemp, 4, 1, F);
	fwrite(pos.x, sizeof(double), n, F);
	fwrite(&ntemp, 4, 1, F);

	fwrite(&ntemp, 4, 1, F);
	fwrite(pos.y, sizeof(double), n, F);
	fwrite(&ntemp, 4, 1, F);

	fwrite(&ntemp, 4, 1, F);
	fwrite(pos.z, sizeof(double), n, F);
	fwrite(&ntemp, 4, 1, F);

	//ntemp = 4*n;
	//fwrite(&ntemp,4,1,F); fwrite(attr1,4,n,F); fwrite(&ntemp,4,1,F);
	//fwrite(&ntemp,4,1,F); fwrite(attr2,4,n,F); fwrite(&ntemp,4,1,F);
	//fwrite(&ntemp,4,1,F); fwrite(attr3,4,n,F); fwrite(&ntemp,4,1,F);
	fclose(F);
}

#include "updateValues.h"

extern struct vec3D pos, momentum;
extern double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE],
		phi[GRID_SIZE][GRID_SIZE][GRID_SIZE],
		delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

void update_density(Shape shape, double a) {
	double dx, dy, dz;
	double tx, ty, tz;

	int x, y, z, i, j, k;
	int n, i_max, j_max, k_max;
	//"reset" the density of the grid

	for (x = 0; x < GRID_SIZE; x++) {
		for (y = 0; y < GRID_SIZE; y++) {
			for (z = 0; z < GRID_SIZE; z++) {
				rho[x][y][z] = 0.0;
			}
		}
	}

	//loop through every particle, index of relevant cell is floor( x/(grid size) )
	for (n = 0; n < N_PARTICLES; n++) {
		i = (int) (pos.x[n]);
		dx = (pos.x[n] - i);
		tx = 1.0 - dx;
		i_max = (i + 1) % GRID_SIZE;

		j = (int) (pos.y[n]);
		dy = pos.y[n] - j;
		ty = 1.0 - dy;
		j_max = (j + 1) % GRID_SIZE;

		k = (int) (pos.z[n]);
		dz = pos.z[n] - k;
		tz = 1.0 - dz;
		k_max = (k + 1) % GRID_SIZE;

		switch (shape) {
		case NGP:
			rho[i][j][k] += M_PARTICLE;
			break;
		case CIC:
			rho[i][j][k] += tx * ty * tz *  M_PARTICLE;
			rho[i_max][j][k] += dx * ty * tz *  M_PARTICLE;
			rho[i][j_max][k] += tx * dy * tz *  M_PARTICLE;
			rho[i][j][k_max] += tx * ty * dz *  M_PARTICLE;
			rho[i][j_max][k_max] += tx * dy * dz *  M_PARTICLE;
			rho[i_max][j_max][k] += dx * dy * tz *  M_PARTICLE;
			rho[i_max][j][k_max] += dx * ty * dz *  M_PARTICLE;
			rho[i_max][j_max][k_max] += dx * dy * dz *  M_PARTICLE;
			break;
		default:
			fprintf(stderr,"Bad shape value (%d) in updateValues::update_density\n",(int) shape);
		}
	}

	for (x = 0; x < GRID_SIZE; x++) {
		for (y = 0; y < GRID_SIZE; y++) {
			for (z = 0; z < GRID_SIZE; z++) {
				delta[x][y][z] = rho[x][y][z] - 1.0;
			}
		}
	}
}

void update_particles(Shape shape, double a, double accelerationSign) {
	const double f_of_a =1.0 / sqrt((OMEGA_M + OMEGA_L * pow(a, 3.0)) / a);
	const double f_of_a_nhalf = 1.0 / sqrt( (OMEGA_M + OMEGA_L * pow((a + DELTA_A / 2.0), 3.0)) / (a + DELTA_A / 2.0));

	//calculate acceleration for each grid
	int n;
	#pragma omp parallel for private(n)
	for (n = 0; n < N_PARTICLES; n++) {
		double accelx, accely, accelz;
		int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max, i_max2, i_min2,
				j_max2, j_min2, k_max2, k_min2;
		double dx, dy, dz;
		double tx, ty, tz;

		i = (int) (pos.x[n]);
		dx = (pos.x[n] - i);
		tx = 1.0 - dx;
		i_max = (i + 1) % GRID_SIZE;
		i_min = (i - 1 + GRID_SIZE) % GRID_SIZE;
		i_min2 = i;
		i_max2 = (i + 2) % GRID_SIZE;

		j = (int) (pos.y[n]);
		dy = pos.y[n] - j;
		ty = 1.0 - dy;
		j_max = (j + 1) % GRID_SIZE;
		j_min = (j - 1 + GRID_SIZE) % GRID_SIZE;
		j_min2 = j;
		j_max2 = (j + 2) % GRID_SIZE;

		k = (int) (pos.z[n]);
		dz = pos.z[n] - k;
		tz = 1.0 - dz;
		k_max = (k + 1) % GRID_SIZE;
		k_min = (k - 1 + GRID_SIZE) % GRID_SIZE;
		k_min2 = k;
		k_max2 = (k + 2) % GRID_SIZE;

		switch (shape) {
		case NGP:
			accelx = accelerationSign * -(phi[i_max][j][k] - phi[i_min][j][k]) / 2.0;
			accely = accelerationSign * -(phi[i][j_max][k] - phi[i][j_min][k]) / 2.0;
			accelz = accelerationSign * -(phi[i][j][k_max] - phi[i][j][k_min]) / 2.0;
			break;
		case CIC:
			accelx = accelerationSign * (-(phi[i_max][j][k] - phi[i_min][j][k]) / 2.0) * tx * ty * tz
					+ (-(phi[i_max2][j][k] - phi[i_min2][j][k]) / 2.0) * dx * ty * tz
					+ (-(phi[i_max][j_max][k] - phi[i_min][j_max][k]) / 2.0) * tx * dy * tz
					+ (-(phi[i_max2][j_max][k] - phi[i_min2][j_max][k]) / 2.0) * dx * dy * tz
					+ (-(phi[i_max][j][k_max] - phi[i_min][j][k_max]) / 2.0) * tx * ty * dz
					+ (-(phi[i_max2][j][k_max] - phi[i_min2][j][k_max]) / 2.0) * dx * ty * dz
					+ (-(phi[i_max][j_max][k_max] - phi[i_min][j_max][k_max]) / 2.0) * tx * dy * dz
					+ (-(phi[i_max2][j_max][k_max] - phi[i_min2][j_max][k_max]) / 2.0) * dx * dy * dz;
			accely = accelerationSign * (-(phi[i][j_max][k] - phi[i][j_min][k]) / 2.0) * tx * ty * tz
					+ (-(phi[i_max][j_max][k] - phi[i_max][j_min][k]) / 2.0) * dx * ty * tz
					+ (-(phi[i][j_max2][k] - phi[i][j_min2][k]) / 2.0) * tx * dy * tz
					+ (-(phi[i_max][j_max2][k] - phi[i_max][j_min2][k]) / 2.0) * dx * dy * tz
					+ (-(phi[i][j_max][k_max] - phi[i][j_min][k_max]) / 2.0) * tx * ty * dz
					+ (-(phi[i_max][j_max][k_max] - phi[i_max][j_min][k_max]) / 2.0) * dx * ty * dz
					+ (-(phi[i][j_max2][k_max] - phi[i][j_min2][k_max]) / 2.0) * tx * dy * dz
					+ (-(phi[i_max][j_max2][k_max] - phi[i_max][j_min2][k_max]) / 2.0) * dx * dy * dz;
			accelz = accelerationSign * (-(phi[i][j][k_max] - phi[i][j][k_min]) / 2.0) * tx * ty * tz
					+ (-(phi[i_max][j][k_max] - phi[i_max][j][k_min]) / 2.0) * dx * ty * tz
					+ (-(phi[i][j_max][k_max] - phi[i][j_max][k_min]) / 2.0) * tx * dy * tz
					+ (-(phi[i_max][j_max][k_max] - phi[i_max][j_max][k_min]) / 2.0) * dx * dy * tz
					+ (-(phi[i][j][k_max2] - phi[i][j][k_min2]) / 2.0) * tx * ty * dz
					+ (-(phi[i_max][j][k_max2] - phi[i_max][j][k_min2]) / 2.0) * dx * ty * dz
					+ (-(phi[i][j_max][k_max2] - phi[i][j_max][k_min2]) / 2.0) * tx * dy * dz
					+ (-(phi[i_max][j_max][k_max2] - phi[i_max][j_max][k_min2]) / 2.0) * dx * dy * dz;

			break;
		     
		default:
			fprintf(stderr, "Bad shape value (%d) in updateValues::update_density\n", (int) shape);
		}
//change momentum first, then position
		momentum.x[n] += f_of_a * accelx * DELTA_A;
		momentum.y[n] += f_of_a * accely * DELTA_A;
		momentum.z[n] += f_of_a * accelz * DELTA_A;

		pos.x[n] += pow((a + DELTA_A / 2.0), -2.0) * f_of_a_nhalf * momentum.x[n] * DELTA_A;
		pos.y[n] += pow((a + DELTA_A / 2.0), -2.0) * f_of_a_nhalf * momentum.y[n] * DELTA_A;
		pos.z[n] += pow((a + DELTA_A / 2.0), -2.0) * f_of_a_nhalf * momentum.z[n] * DELTA_A;

		if (pos.x[n] > L_BOX)
			pos.x[n] -= L_BOX;
		if (pos.y[n] > L_BOX)
			pos.y[n] -= L_BOX;
		if (pos.z[n] > L_BOX)
			pos.z[n] -= L_BOX;

		if (pos.x[n] < 0.)
			pos.x[n] += L_BOX;
		if (pos.y[n] < 0.)
			pos.y[n] += L_BOX;
		if (pos.z[n] < 0.)
			pos.z[n] += L_BOX;
	}
}

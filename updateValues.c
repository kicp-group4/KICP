#include "updateValues.h"
#include <sys/time.h>

extern struct vec3D pos, momentum;
extern double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE],
		phi[GRID_SIZE][GRID_SIZE][GRID_SIZE],
		delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

void update_density(double a) {
	double dx, dy, dz;
	double tx, ty, tz;

	int i, j, k;
	int n, i_max, j_max, k_max;
	//"reset" the density of the grid

	double *p_rho = &rho[0][0][0];
	#pragma omp parallel for private(n)
	for(n=0;n<GRID_SIZE*GRID_SIZE*GRID_SIZE;n+=2){
		p_rho[n] = 0.0;
		p_rho[n+1] = 0.0;
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

#ifdef SF_NGP
		rho[i][j][k] += M_PARTICLE;
#elif defined SF_CIC
		rho[i][j][k] += tx * ty * tz;
		rho[i_max][j][k] += dx * ty * tz;
		rho[i][j_max][k] += tx * dy * tz;
		rho[i][j][k_max] += tx * ty * dz;
		rho[i][j_max][k_max] += tx * dy * dz;
		rho[i_max][j_max][k] += dx * dy * tz;
		rho[i_max][j][k_max] += dx * ty * dz;
		rho[i_max][j_max][k_max] += dx * dy * dz;
#else
#		error "Must specify a shape function"
#endif
	}

	double *p_delta = &delta[0][0][0];
	p_rho = &rho[0][0][0];
	#pragma omp parallel for private(n)
	for(n=0;n<GRID_SIZE*GRID_SIZE*GRID_SIZE;n+=2){
		p_rho[n] *= M_PARTICLE;
		p_rho[n+1] *= M_PARTICLE;
		p_delta[n] = p_rho[n] - 1.0;
		p_delta[n+1] = p_rho[n+1] - 1.0;
	}
}

void update_particles(const double a) {
	const register double f_of_a =1.0 / sqrt((OMEGA_M + OMEGA_L * pow(a, 3.0)) / a);
	const register double f_of_a_nhalf = 1.0 / sqrt( (OMEGA_M + OMEGA_L * pow((a + DELTA_A / 2.0), 3.0)) / (a + DELTA_A / 2.0));
	const register double momentumScale = f_of_a * DELTA_A;
	const register double positionScale = f_of_a_nhalf * DELTA_A * pow((a + DELTA_A / 2.0), -2.0);

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

#ifdef SF_NGP
		accelx = ST_GLASS_ACCEL_SIGN -(phi[i_max][j][k] - phi[i_min][j][k]) / 2.0;
		accely = ST_GLASS_ACCEL_SIGN -(phi[i][j_max][k] - phi[i][j_min][k]) / 2.0;
		accelz = ST_GLASS_ACCEL_SIGN -(phi[i][j][k_max] - phi[i][j][k_min]) / 2.0;
#elif defined SF_CIC
		accelx = (-(phi[i_max][j][k] - phi[i_min][j][k]) / 2.0) * tx * ty * tz
				+ (-(phi[i_max2][j][k] - phi[i_min2][j][k]) / 2.0) * dx * ty * tz
				+ (-(phi[i_max][j_max][k] - phi[i_min][j_max][k]) / 2.0) * tx * dy * tz
				+ (-(phi[i_max2][j_max][k] - phi[i_min2][j_max][k]) / 2.0) * dx * dy * tz
				+ (-(phi[i_max][j][k_max] - phi[i_min][j][k_max]) / 2.0) * tx * ty * dz
				+ (-(phi[i_max2][j][k_max] - phi[i_min2][j][k_max]) / 2.0) * dx * ty * dz
				+ (-(phi[i_max][j_max][k_max] - phi[i_min][j_max][k_max]) / 2.0) * tx * dy * dz
				+ (-(phi[i_max2][j_max][k_max] - phi[i_min2][j_max][k_max]) / 2.0) * dx * dy * dz;
		accely = (-(phi[i][j_max][k] - phi[i][j_min][k]) / 2.0) * tx * ty * tz
				+ (-(phi[i_max][j_max][k] - phi[i_max][j_min][k]) / 2.0) * dx * ty * tz
				+ (-(phi[i][j_max2][k] - phi[i][j_min2][k]) / 2.0) * tx * dy * tz
				+ (-(phi[i_max][j_max2][k] - phi[i_max][j_min2][k]) / 2.0) * dx * dy * tz
				+ (-(phi[i][j_max][k_max] - phi[i][j_min][k_max]) / 2.0) * tx * ty * dz
				+ (-(phi[i_max][j_max][k_max] - phi[i_max][j_min][k_max]) / 2.0) * dx * ty * dz
				+ (-(phi[i][j_max2][k_max] - phi[i][j_min2][k_max]) / 2.0) * tx * dy * dz
				+ (-(phi[i_max][j_max2][k_max] - phi[i_max][j_min2][k_max]) / 2.0) * dx * dy * dz;
		accelz = (-(phi[i][j][k_max] - phi[i][j][k_min]) / 2.0) * tx * ty * tz
				+ (-(phi[i_max][j][k_max] - phi[i_max][j][k_min]) / 2.0) * dx * ty * tz
				+ (-(phi[i][j_max][k_max] - phi[i][j_max][k_min]) / 2.0) * tx * dy * tz
				+ (-(phi[i_max][j_max][k_max] - phi[i_max][j_max][k_min]) / 2.0) * dx * dy * tz
				+ (-(phi[i][j][k_max2] - phi[i][j][k_min2]) / 2.0) * tx * ty * dz
				+ (-(phi[i_max][j][k_max2] - phi[i_max][j][k_min2]) / 2.0) * dx * ty * dz
				+ (-(phi[i][j_max][k_max2] - phi[i][j_max][k_min2]) / 2.0) * tx * dy * dz
				+ (-(phi[i_max][j_max][k_max2] - phi[i_max][j_max][k_min2]) / 2.0) * dx * dy * dz;
#else
#	error "Must specify a shape function"
#endif

#ifdef ST_GLASS
	#define ST_GLASS_ACCEL_SIGN -
#else
	#define ST_GLASS_ACCEL_SIGN
#endif

//change momentum first, then position
		momentum.x[n] += ST_GLASS_ACCEL_SIGN accelx * momentumScale;
		momentum.y[n] += ST_GLASS_ACCEL_SIGN accely * momentumScale;
		momentum.z[n] += ST_GLASS_ACCEL_SIGN accelz * momentumScale;

#undef ST_GLASS_ACCEL_SIGN

		pos.x[n] += momentum.x[n] * positionScale;
		pos.y[n] += momentum.y[n] * positionScale;
		pos.z[n] += momentum.z[n] * positionScale;

		/* Enforce periodic boundary conditions */
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

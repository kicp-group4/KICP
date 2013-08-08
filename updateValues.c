#include "updateValues.h"

extern struct vec3D pos, momentum;
extern double rho[GRID_SIZE][GRID_SIZE][GRID_SIZE],
	phi[GRID_SIZE][GRID_SIZE][GRID_SIZE],
	delta[GRID_SIZE][GRID_SIZE][GRID_SIZE];

void update_density(Shape shape, double a) {
	double dx, dy, dz;
	double tx, ty, tz;
	double a_cubed = a*a*a;
	double scale = DELTA_X / a;
	int x, y, z, i, j, k;
	int n,i_max,j_max,k_max;
	double R_0_3 = R_0*R_0*R_0;
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
	
	  i = floor((pos.x[n]/scale) - 0.5);
		dx = pos.x[n]/scale - i ;
		tx = 1.0 - dx;
		if (i < 0){
		  i = GRID_SIZE - 1;
		  dx = (pos.x[n]/scale);
		  tx = 1.0 - dx;
		}
		i_max = i + 1;
		if (i_max == GRID_SIZE) i_max = 0;
		
		j = floor((pos.y[n] / scale) - 0.5);
		dy = pos.y[n]/scale - j ;
		ty = 1.0 - dy;		
		if (j < 0){
		  j = GRID_SIZE - 1;
		  dy = (pos.y[n]/scale);
		  ty = 1.0 - dy;
		}
		j_max = j + 1;
		if (j_max == GRID_SIZE) j_max = 0;
		
		k = floor((pos.z[n] / scale) - 0.5);
		dz = pos.z[n]/scale - k;
		tz = 1.0 - dz;
		if (k < 0){
		  k =  GRID_SIZE - 1;
		  dz = (pos.z[n]/scale);
		  tz = 1.0 - dz;
		}
		k_max = k + 1;
		if (k_max == GRID_SIZE) k_max = 0;

		/* printf("%g\t%g\t%g\n",pos.x[n],pos.y[n],pos.z[n]); */
		/* printf("%d\t%d\t%d\n",i,j,k); */
		

		switch(shape) {
		case NGP:
			rho[i][j][k] += a_cubed*(M_PARTICLE * pow((1.0 / scale), 3.0)) / (RHO_0);
			break;
		case CIC:
			rho[i][j][k] += a_cubed*(M_PARTICLE/R_0_3 * tx * ty * tz) / RHO_0;
			rho[i_max][j][k] += a_cubed*(M_PARTICLE/R_0_3 * dx * ty * tz) / RHO_0;
			rho[i][j_max][k] += a_cubed*(M_PARTICLE/R_0_3 * tx * dy * tz) / RHO_0;
			rho[i][j][k_max] += a_cubed*(M_PARTICLE/R_0_3 * tx * ty * dz) / RHO_0;
			rho[i][j_max][k_max] += a_cubed*(M_PARTICLE/R_0_3 * tx * dy * dz) / RHO_0;
			rho[i_max][j_max][k] += a_cubed*(M_PARTICLE/R_0_3 * dx * dy * tz) / RHO_0;
			rho[i_max][j][k_max] += a_cubed*(M_PARTICLE/R_0_3 * dx * ty * dz) / RHO_0;
			rho[i_max][j_max][k_max] += a_cubed*(M_PARTICLE/R_0_3 * dx * dy * dz) / RHO_0;

			break;
		default:
			fprintf(stderr,"Bad shape value (%d) in updateValues::update_density\n",(int)shape);
		}
	}

	for(x=0;x<GRID_SIZE;x++) {
		for(y=0;y<GRID_SIZE;y++) {
			for(z=0;z<GRID_SIZE;z++) {
				delta[x][y][z] = rho[x][y][z] - 1.0;
			}
		}
	}
}

void update_particles(Shape shape, double a) {
	double f_of_a, f_of_a_nhalf;
	double accelx, accely, accelz;
	int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;
	double dx, dy, dz;
	double tx, ty, tz;
	double scale = DELTA_X / a;
	double Lb = L_BOX / (R_0 * a);


	f_of_a = 1.0 / sqrt((OMEGA_M + OMEGA_L * pow(a, 3.0)) / a);
	f_of_a_nhalf = 1.0/ sqrt((OMEGA_M + OMEGA_L * pow((a + DELTA_A / 2.0), 3.0))/ (a + DELTA_A / 2.0));

	//calculate acceleration for each grid
	int n;
	for (n = 0; n < N_PARTICLES; n++) {
		i = floor((pos.x[n]/scale) - 0.5);
		dx = pos.x[n]/scale - i ;
		tx = 1.0 - dx;
		if (i < 0){
		  i = GRID_SIZE - 1;
		  dx = (pos.x[n]/scale);
		  tx = 1.0 - dx;
		}
		i_min = i - 1;
		i_max = i + 1;
		if (i_min < 0) i_min = GRID_SIZE - 1;
		if (i_max == GRID_SIZE) i_max = 0;
		
		j = floor((pos.y[n] / scale) - 0.5);
		dy = pos.y[n]/scale - j ;
		ty = 1.0 - dy;		
		if (j < 0){
		  j = GRID_SIZE - 1;
		  dy = (pos.y[n]/scale);
		  ty = 1.0 - dy;
		}
		j_min = j - 1;
		j_max = j + 1;
		if (j_min < 0) j_min = GRID_SIZE - 1;
		if (j_max == GRID_SIZE) j_max = 0;
		
		k = floor((pos.z[n] / scale) - 0.5);
		dz = pos.z[n]/scale - k;
		tz = 1.0 - dz;
		if (k < 0){
		  k =  GRID_SIZE - 1;
		  dz = (pos.z[n]/scale);
		  tz = 1.0 - dz;
		}
		k_min = k - 1;
		k_max = k + 1;
		if (k_min < 0) k_min = GRID_SIZE - 1;
		if (k_max == GRID_SIZE) k_max = 0;

		switch(shape) {
		case NGP:
			accelx = -(phi[i_max][j][k] - phi[i_min][j][k]) / 2.0;
			accely = -(phi[i][j_max][k] - phi[i][j_min][k]) / 2.0;
			accelz = -(phi[i][j][k_max] - phi[i][j][k_min]) / 2.0;
			break;
		case CIC:
			accelx = (-(phi[i_max][j][k] - phi[i_min][j][k]) / 2.0) * tx * ty * tz
					+ (-(phi[i_max + 1][j][k] - phi[i_min + 1][j][k]) / 2.0) * dx * ty * tz
					+ (-(phi[i_max][j + 1][k] - phi[i_min][j + 1][k]) / 2.0) * tx * dy * tz
					+ (-(phi[i_max + 1][j + 1][k] - phi[i_min + 1][j + 1][k]) / 2.0) * dx * dy * tz
					+ (-(phi[i_max][j][k + 1] - phi[i_min][j][k + 1]) / 2.0) * tx * ty * dz
					+ (-(phi[i_max + 1][j][k + 1] - phi[i_min + 1][j][k + 1]) / 2.0) * dx * ty * dz
					+ (-(phi[i_max][j + 1][k + 1] - phi[i_min][j + 1][k + 1]) / 2.0) * tx * dy * dz
					+ (-(phi[i_max + 1][j + 1][k + 1] - phi[i_min + 1][j + 1][k + 1]) / 2.0) * dx * dy * dz;
			accely = (-(phi[i][j_max][k] - phi[i][j_min][k]) / 2.0) * tx * ty * tz
					+ (-(phi[i + 1][j_max][k] - phi[i + 1][j_min][k]) / 2.0) * dx * ty * tz
					+ (-(phi[i][j_max + 1][k] - phi[i][j_min + 1][k]) / 2.0) * tx * dy * tz
					+ (-(phi[i + 1][j_max + 1][k] - phi[i + 1][j_min + 1][k]) / 2.0) * dx * dy * tz
					+ (-(phi[i][j_max][k + 1] - phi[i][j_min][k + 1]) / 2.0) * tx * ty * dz
					+ (-(phi[i + 1][j_max][k + 1] - phi[i + 1][j_min][k + 1]) / 2.0) * dx * ty * dz
					+ (-(phi[i][j_max + 1][k + 1] - phi[i][j_min + 1][k + 1]) / 2.0) * tx * dy * dz
					+ (-(phi[i + 1][j_max + 1][k + 1] - phi[i + 1][j_min + 1][k + 1]) / 2.0) * dx * dy * dz;
			accelz = (-(phi[i][j][k_max] - phi[i][j][k_min]) / 2.0) * tx * ty * tz
					+ (-(phi[i + 1][j][k_max] - phi[i + 1][j][k_min]) / 2.0) * dx * ty * tz
					+ (-(phi[i][j + 1][k_max] - phi[i][j + 1][k_min]) / 2.0) * tx * dy * tz
					+ (-(phi[i + 1][j + 1][k_max] - phi[i + 1][j + 1][k_min]) / 2.0) * dx * dy * tz
					+ (-(phi[i][j][k_max + 1] - phi[i][j][k_min + 1]) / 2.0) * tx * ty * dz
					+ (-(phi[i + 1][j][k_max + 1] - phi[i + 1][j][k_min + 1]) / 2.0) * dx * ty * dz
					+ (-(phi[i][j + 1][k_max + 1] - phi[i][j + 1][k_min + 1]) / 2.0) * tx * dy * dz
					+ (-(phi[i + 1][j + 1][k_max + 1] - phi[i + 1][j + 1][k_min + 1]) / 2.0) * dx * dy * dz;
			break;
		default:
			fprintf(stderr,"Bad shape value (%d) in updateValues::update_density\n",(int)shape);
		}
//change momentum first, then position
		momentum.x[n] += f_of_a * accelx * DELTA_A;
		momentum.y[n] += f_of_a * accely * DELTA_A;
		momentum.z[n] += f_of_a * accelz * DELTA_A;


		pos.x[n] += pow((a + DELTA_A / 2.0), -2.0) * f_of_a_nhalf * momentum.x[n] * DELTA_A;
		pos.y[n] += pow((a + DELTA_A / 2.0), -2.0) * f_of_a_nhalf * momentum.y[n] * DELTA_A;
		pos.z[n] += pow((a + DELTA_A / 2.0), -2.0) * f_of_a_nhalf * momentum.z[n] * DELTA_A;

		if(pos.x[n] > Lb) pos.x[n] -= Lb;
		if(pos.y[n] > Lb) pos.y[n] -= Lb;
		if(pos.z[n] > Lb) pos.z[n] -= Lb;
		  
		if(pos.x[n] < 0.) pos.x[n] += Lb;
		if(pos.y[n] < 0.) pos.y[n] += Lb;
		if(pos.z[n] < 0.) pos.z[n] += Lb;
	}
}

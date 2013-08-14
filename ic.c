#include "ic.h"

extern struct vec3D pos, momentum;

void ic(double a) {
	int i, j, m, Np;
	double q, D_aini, D_dot, k;
	/*!
	 * 
	 * 	1 - D
	 * 
	 * 	C O L L A P S E
	 * 
	 * 	P L A N E  W A V E 
	 * 
	 */
	// POSITION + VELOCITY 
	/* NOTE: D_+(a)=a/a_0, but we assume a_0=1 throughout */
	Np = N_P_1D;
	k = 2 * M_PI / L_BOX;
	D_aini = a;

	double a_cross = 10.0 * A_INITIAL;
	double amp = 1.0 / (a_cross * k);

	double a_temp = A_INITIAL - DELTA_A / 2.0;

	D_dot = sqrt(cosmology->OmegaM + cosmology->OmegaL * pow(a_temp, 3.0)) / sqrt(a_temp);

	for (i = 0; i < Np; i++) {
	  q = i*(GRID_SIZE/Np);
		for (j = 0; j < Np; j++) {
			for (m = 0; m < Np; m++) {
				/* NOTE: 2*M_PI*q/L_BOX = k*q = 2*M_PI*i/Np */
				pos.x[i + Np * (j + Np * m)] = (q + D_aini * amp * sin(2 * M_PI * i / Np));

				if (pos.x[i + Np * (j + Np * m)] >= GRID_SIZE) {
					pos.x[i + Np * (j + Np * m)] = pos.x[i + Np * (j + Np * m)] - L_BOX;
				}
				if (pos.x[i + Np * (j + Np * m)] < 0) {
					pos.x[i + Np * (j + Np * m)] = pos.x[i + Np * (j + Np * m)] + L_BOX;
				}

				momentum.x[i + Np * (j + Np * m)] = a_temp * (a_temp * D_dot * amp * sin(2 * M_PI * i / Np));
				pos.y[i + Np * (j + Np * m)] = j*(GRID_SIZE/Np);
				momentum.y[i + Np * (j + Np * m)] = 0;
				pos.z[i + Np * (j + Np * m)] = m*(GRID_SIZE/Np);
				momentum.z[i + Np * (j + Np * m)] = 0;
				
			}
		}
	       
	}

}

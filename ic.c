#include "ic.h"

extern struct vec3D pos, momentum;

void ic(double a) {
	int i, j, m, Np;
	double q, D_aini, D_dot, k, L_NP;
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
	Np = 32;
	L_NP = L_BOX / (Np-1);
	k = 2 * PI / L_BOX;
	D_aini = 1;

	double a_cross = 10.0 * A_INITIAL;
	double amp = -1.0 / (a_cross * k / A_INITIAL);

	double a_temp = A_INITIAL - DELTA_A / 2.0;
	D_dot = ((H_0 * 3.241e-20) * sqrt(OMEGA_M + OMEGA_L * pow(a_temp, 3.0)) / sqrt(a_temp)) / A_INITIAL;

	for (i = 0; i < Np; i++) {
	  q = i * L_NP;
		for (j = 0; j < Np; j++) {
			for (m = 0; m < Np; m++) {
				pos.x[i + Np * (j + Np * m)] = (q + D_aini * amp * sin(k * q)) / R_0 / a;
				momentum.x[i + Np * (j + Np * m)] = a*(A_INITIAL * D_dot * amp * sin(k * q)) / (V_0);
				pos.y[i + Np * (j + Np * m)] = (j * L_NP) / R_0 / a;
				momentum.y[i + Np * (j + Np * m)] = 0;
				pos.z[i + Np * (j + Np * m)] = (m * L_NP) / R_0 / a;
				momentum.z[i + Np * (j + Np * m)] = 0;
			}
		}
	}
	
	
// 	pos.x[0] = 8.0 * L_BOX / GRID_SIZE / R_0 / a;
// 	pos.y[0] = 0.0;
// 	pos.z[0] = 0.0;
// 	pos.x[1] = 24.0 * L_BOX / GRID_SIZE / R_0 / a;
// 	pos.y[1] = 0.0;
// 	pos.z[1] = 0.0;
// 	momentum.x[0] = 0.0;
// 	momentum.y[0] = 0.; //0.5*sqrt(M_PARTICLE*G/(16.0*L_NP))/(M_PARTICLE*V_0);
// 	momentum.z[0] = 0.; //0.5*sqrt(M_PARTICLE*G/(16.0*L_NP))/(M_PARTICLE*V_0);
// 	momentum.x[1] = 0.0;
// 	momentum.y[1] = 0.; //-0.5*sqrt(M_PARTICLE*G/(16.0*L_NP))/(M_PARTICLE*V_0);
// 	momentum.z[1] = 0.; //-0.5*sqrt(M_PARTICLE*G/(16.0*L_NP))/(M_PARTICLE*V_0);
}

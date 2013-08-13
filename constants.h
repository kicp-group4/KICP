#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define PI 3.14159

/* ***************************************************
 * Make sure each of these is a power of 2. If they
 * are not, you are going to have a very bad time.
 * */
#define GRID_SIZE 256
#define L_BOX 256.0		/* This is in MegaParsecs*/
#define N_P_1D 64
/****************************************************/

#define N_PARTICLES ((N_P_1D)*(N_P_1D)*(N_P_1D))
#define G 6.67384e-11
#define M_PARTICLE ((GRID_SIZE*GRID_SIZE*GRID_SIZE)/((float)N_PARTICLES))

#ifdef ST_COSMO_ZELDOVICH
	#define OMEGA_M 0.275
	#define OMEGA_L 0.725
	#define DELTA_A 0.01
	#define A_INITIAL (1/(1.0+Z_INI))
#else
	#define OMEGA_M 1.0
	#define OMEGA_L 0.0
	#define DELTA_A 0.01
	#define A_INITIAL 0.2
#endif

//dimensionless scaling values

#define R_0 (L_BOX/GRID_SIZE)
#define V_0 (R_0/T_0)
#define RHO_0 (3*OMEGA_M/(8*PI*G*T_0*T_0))
#define PHI_0 (V_0*V_0)

#endif

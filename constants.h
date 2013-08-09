#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define PI 3.14159

#define GRID_SIZE 64
#define L_BOX 64.0
#define N_P_1D 32
#define N_PARTICLES (N_P_1D)*(N_P_1D)*(N_P_1D)
#define H_0 (72*3.241e-20)
#define G 6.67384e-11
#define OMEGA_M 1.0
#define OMEGA_L 0.0
#define M_PARTICLE ((L_BOX*L_BOX*L_BOX)/((float)N_PARTICLES))

#define DELTA_T 0.01
#define DELTA_A 0.02
#define A_INITIAL 0.2
#define DELTA_X 1.0

//dimensionless scaling values

#define R_0 (L_BOX/GRID_SIZE)
#define T_0 (1/H_0)
#define V_0 (R_0/T_0)
#define RHO_0 (3*OMEGA_M/(8*PI*G*T_0*T_0))
#define PHI_0 (V_0*V_0)

#endif

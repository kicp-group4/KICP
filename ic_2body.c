//#include "ic.h"
//
//extern struct vec3D pos, momentum;
//
//void ic_2body(double a){
//
//	FILE *ic;
//	int i, j, m, Np;
//	double q, D_aini, D_dot, k, L_NP;
//	char ic_filename[] = "ic_2body.dat";
//	ic = fopen(ic_filename,"w");
//
//	/*
//	 * HERE'S THE CODE FOR THE HEADER IF YOU WANT TO INCLUDE IT LATER
//	 *
//	// HEADER
//	fprintf(output,"%g\n",L_BOX);
//	fprintf(output,"%d\n",GRID_SIZE);
//	fprintf(output,"%d\n",N_PARTICLES);
//	fprintf(output,"%g\n",H_0);
//	fprintf(output,"%g\n",OMEGA_M);
//	fprintf(output,"%g\n",OMEGA_L);
//	fprintf(output,"%g\n",M_PARTICLE);
//	*/
//	/*!
//	 *
//	 * 	2-BODY ORBIT -> MAKE SURE TO CHANGE N_PARTICLES TO 2, delta_a needs to be small b/c momentum is not offset
//	 *
//	 *
//	 */
//	// POSITION + VELOCITY
//
//	pos.x[0] = (8.0 * L_BOX/GRID_SIZE)/R_0
//	momentum.x[0] = 0.0;
//	pos.y[0] = 0.0;
//	momentum.y[0] = sqrt(M_PARTICLE*G/(8.0 * L_BOX/GRID_SIZE))/(M_PARTICLE*V_0);
//	pos.z[0] = 0.0;
//	momentum.z[0] = 0.0;
//	fprintf(ic, "%g \t %g \t %g", pos.x[0], pos.y[0], pos.z[0]);
//	fprintf(ic," \t %g \t %g \t %g \n", momentum.x[0],momentum.y[0],momentum.z[0]);
//
//	pos.x[1] = (24.0 * L_BOX/GRID_SIZE)/R_0
//	momentum.x[1] = 0.0;
//	pos.y[1] = 0.0;
//	momentum.y[1] = -1.0*sqrt(M_PARTICLE*G/(8.0 * L_BOX/GRID_SIZE))/(M_PARTICLE*V_0);
//	pos.z[1] = 0.0;
//	momentum.z[1] = 0.0;
//	fprintf(ic, "%g \t %g \t %g", pos.x[1], pos.y[1], pos.z[1]);
//	fprintf(ic," \t %g \t %g \t %g \n", momentum.x[1],momentum.y[1],momentum.z[1]);
//
//	fclose(ic);
//}

#ifndef POISSONSOLVER_H_
#define POISSONSOLVER_H_
#include "constants.h"
#include <fftw3.h>
#include <math.h>

void poissonSolver(double a);
void poissonSolver_init(void);
void poissonSolver_cleanup(void);


#endif /* POISSONSOLVER_H_ */

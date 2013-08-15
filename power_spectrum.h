#ifndef _POWER_SPECTRUM_H_
#define _POWER_SPECTRUM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "vec3D.h"
#include "cosmology.h"
#define NUM_BINS GRID_SIZE

void power_spectrum(double);
void power_spectrum_init();
void power_spectrum_cleanup();

#endif

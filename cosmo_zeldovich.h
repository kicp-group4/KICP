#ifndef _COSMO_ZELDOVICH_H_
#define _COSMO_ZELDOVICH_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "vec3D.h"

void cosmo_zeldovich(double a);
void cosmo_zeldovich_init(int n_ks);
void cosmo_zeldovich_cleanup();

#endif
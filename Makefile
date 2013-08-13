CC=icc
CFLAGS = -c -O3 -Wall -openmp -DSF_CIC -DST_COSMO_ZELDOVICH -DZ_INI=100
LDFLAGS= -openmp -lm -lfftw3_omp -lfftw3 -lgsl -lgslcblas
SOURCES= cosmo_zeldovich.c cosmology.c ic.c ic_glass.c output.c pm.c poissonSolver.c power_spectrum.c updateValues.c 
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=pm

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) -o $@ $(CFLAGS) $<

clean:
	rm -f *.o *.dat pm test.out test.err gmon.out

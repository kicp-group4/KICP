CC=icc
CFLAGS = -c -O1 -Wall -openmp -p -DSF_CIC -DST_ZELDOVICH
LDFLAGS= -openmp -p -lm -lfftw3 -lgsl -lgslcblas
SOURCES=pm.c ic.c output.c poissonSolver.c updateValues.c ic_glass.c power_spectrum.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=pm

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) -o $@ $(CFLAGS) $<

clean:
	rm -f *.o *.dat

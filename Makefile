
objects = pm.o ic.o output.o poissonSolver.o updateValues.o

PM_G4 : $(objects)
	gcc  $(objects) -L/usr/lib -lfftw3 -lm -o PM_G4

pm.o : pm.c constants.h vec3D.h ic.h shape.h updateValues.h poissonSolver.h
	gcc -c pm.c

ic.o : ic.c ic.h constants.h vec3D.h
	gcc -c ic.c

output.o : output.c output.h vec3D.h
	gcc -c output.c

poissonSolver.o : poissonSolver.c poissonSolver.h constants.h
	gcc -c poissonSolver.c

updateValues.o : updateValues.c updateValues.h vec3D.h shape.h
	gcc -c updateValues.c

clean:
	rm -f $(objects) $(objetcs_ps) *.dat
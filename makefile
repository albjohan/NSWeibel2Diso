CC=mpicc
CFILES= main.c newton.c
HFILES = newton.h


weib_roots : $(CFILES) $(HFILES)
	$(CC) -march=native -O2 main.c newton.c -o weib_roots -lm -lcerf -lmpi

clean :
	rm weib_roots

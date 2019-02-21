CCOMP = gcc
CPLUSCOMP = g++
MEX=$(MATLAB)/bin/mex

all: flow farey cutweird cutexp
flow: Pairing.c Sparse.c flow.c timer.c types.h timer.h
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) -largeArrayDims Pairing.c Sparse.c flow.c timer.c
farey: Farey.c
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) -largeArrayDims Farey.c 
cutweird: cutweird.c
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) -largeArrayDims cutweird.c
cutexp: cutexp.c
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) -largeArrayDims cutexp.c

clean:
	rm -f *~ 
	rm -f \#*

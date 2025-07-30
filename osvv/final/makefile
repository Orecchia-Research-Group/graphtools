CCOMP=gcc
CPLUSCOMP=g++
MEX=mex

CXXFLAGS=-std=c++20 -fPIC -Wall -Wextra -O3
MEXFLAGS=-O


all: Pairing.mexa64 cutexp.mexa64 cutweird.mexa64 Farey.mexa64 loadMetisGraph.mexa64 loadMtxGraph.mexa64 readPtn.mexa64

debug: set-debug all
set-debug:
	MEXFLAGS=-g -DDEBUG
Pairing.mexa64: Pairing.o Sparse.o flow.o timer.o dynamictree.o
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^
cutweird.mexa64: cutweird.o farey.o
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^
Farey.mexa64: Farey.o farey.o
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^
%.mexa64: %.c
	$(MEX) CC=$(CCOMP) LD=$(CPLUSCOMP) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^
%.mexa64: %.cpp
	$(MEX) CPP=$(CPLUSCOMP) LD=$(CPLUSCOMP) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

%.o: %.c
	$(MEX) CC=$(CCOMP) LD=$(CCOMP) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) -c $^
%.o: %.cpp
	$(MEX) CXX=$(CPLUSCOMP) LD=$(CPLUSCOMP) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) -c $^

clean:
	rm -rf *~ *.mexa64 *.o
	rm -f \#*

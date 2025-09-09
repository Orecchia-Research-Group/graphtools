CC     ?= gcc
CXX    ?= g++
MEX    ?= mex

CXXFLAGS = -std=c++20 -fPIC -Wall -Wextra -O3
MEXFLAGS = -O

# Detect OS and choose mex extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    MEXEXT = mexa64
    OBJEXT = o
else ifeq ($(UNAME_S),Darwin)
    MEXEXT = mexmaci64
    OBJEXT = o
else
    # Windows: default mex extension
    MEXEXT = mexw64
	MEX = mex.bat
    # Use mex.bat -setup to detect the active compiler
    MEXINFO := $(shell $(MEX) -setup 2>&1 | head -n 1)
    ifneq (,$(findstring Microsoft,$(MEXINFO)))
        OBJEXT = obj   # MSVC
        CC = cl
        CXX = cl
    else ifneq (,$(findstring MinGW,$(MEXINFO)))
        OBJEXT = o     # MinGW
        CC = gcc
        CXX = g++
    else ifneq (,$(findstring clang,$(MEXINFO)))
        OBJEXT = o     # Clang
        CC = clang
        CXX = clang++
    else
        # fallback if detection fails
        OBJEXT = obj
    endif
endif

TARGETS = \
    Pairing.$(MEXEXT) \
    cutexp.$(MEXEXT) \
    cutweird.$(MEXEXT) \
    Farey.$(MEXEXT) \
    loadMetisGraph.$(MEXEXT) \
    loadMtxGraph.$(MEXEXT) \
    readPtn.$(MEXEXT)

all: $(TARGETS)

debug: MEXFLAGS=-g -DDEBUG
debug: clean all

Pairing.$(MEXEXT): Pairing.o Sparse.o flow.o timer.o dynamictree.o
	$(MEX) CC=$(CC) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

cutweird.$(MEXEXT): cutweird.o farey.o
	$(MEX) CC=$(CC) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

Farey.$(MEXEXT): Farey.o farey.o
	$(MEX) CC=$(CC) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

%.${MEXEXT}: %.c
	$(MEX) CC=$(CC) LD=$(CC) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

%.${MEXEXT}: %.cpp
	$(MEX) CXX=$(CXX) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

%.o: %.c
	$(MEX) CC=$(CC) LD=$(CC) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) -c $^

%.o: %.cpp
	$(MEX) CXX=$(CXX) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) -c $^

clean:
	$(RM) -f *~ *.o *.$(MEXEXT) \#*

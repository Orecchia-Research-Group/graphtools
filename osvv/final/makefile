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
    $(info $(MEXIFNO))
    ifneq (,$(findstring Microsoft,$(MEXINFO)))
        OBJEXT = obj   # MSVC
        CC = cl
        CXX = cl
    else ifneq (,$(findstring MinGW,$(MEXINFO)))
        OBJEXT = obj     # MinGW
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

Pairing.$(MEXEXT): Pairing.$(OBJEXT) Sparse.$(OBJEXT)  flow.$(OBJEXT) timer.$(OBJEXT) dynamictree.$(OBJEXT)
	$(MEX) CC=$(CC) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

cutweird.$(MEXEXT): cutweird.$(OBJEXT) farey.$(OBJEXT)
	$(MEX) CC=$(CC) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

Farey.$(MEXEXT): Farey.$(OBJEXT) farey.$(OBJEXT)
	$(MEX) CC=$(CC) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

%.${MEXEXT}: %.c
	$(MEX) CC=$(CC) LD=$(CC) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

%.${MEXEXT}: %.cpp
	$(MEX) CXX=$(CXX) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) $^

%.$(OBJEXT): %.c
	$(MEX) CC=$(CC) LD=$(CC) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) -c $^

%.$(OBJEXT): %.cpp
	$(MEX) CXX=$(CXX) LD=$(CXX) CXXFLAGS="$(CXXFLAGS)" $(MEXFLAGS) -c $^

clean:
	$(RM) -f *~ *.$(OBJEXT) *.$(MEXEXT) \#*

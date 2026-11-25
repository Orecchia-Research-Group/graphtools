# Standalone tests for flow.c / dynamictree.c (no MATLAB needed).
#   make          build test_matching
#   make test     correctness sweep: dinic and dynamic on random and grid graphs
#   make bench    timing comparison on larger random and grid graphs
#   make debug    check that a -DDEBUG build of flow.c still links
SRC    = ..
CC    ?= gcc
CXX   ?= g++
CFLAGS = -O2 -std=gnu11 -I$(SRC)

test_matching: test_matching.o flow.o dynamictree.o timer.o
	$(CXX) $^ -o $@

test_matching.o: test_matching.c $(SRC)/flow.h
	$(CC) $(CFLAGS) -c $< -o $@
flow.o: $(SRC)/flow.c $(SRC)/flow.h $(SRC)/dynamictree.h $(SRC)/timer.h
	$(CC) $(CFLAGS) -w -c $< -o $@
dynamictree.o: $(SRC)/dynamictree.c $(SRC)/dynamictree.h
	$(CC) $(CFLAGS) -w -c $< -o $@
timer.o: $(SRC)/timer.cpp $(SRC)/timer.h
	$(CXX) -O2 -I$(SRC) -c $< -o $@

test: test_matching
	@fail=0; for seed in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do \
	  for algo in 0 1; do \
	    for g in "random 6" "random 50" "random 2000" "grid 4" "grid 30"; do \
	      ./test_matching $$g $$seed $$algo > /dev/null 2>&1 || { fail=$$((fail+1)); echo "FAIL: $$g seed=$$seed algo=$$algo"; }; \
	    done; done; done; \
	echo "$$fail failures out of 200 runs"; [ $$fail -eq 0 ]

bench: test_matching
	@for g in "random 200000" "random 1000000" "grid 300" "grid 1000"; do \
	  for algo in -1 0 1; do ./test_matching $$g 1 $$algo 2>/dev/null; done; done

debug: test_matching.c
	$(CC) $(CFLAGS) -w -DDEBUG -c $(SRC)/flow.c -o flow_debug.o
	$(CC) $(CFLAGS) -w -c $(SRC)/dynamictree.c -o dynamictree_debug.o
	$(CC) $(CFLAGS) -c test_matching.c -o test_matching_debug.o
	$(CXX) -O2 -I$(SRC) -c $(SRC)/timer.cpp -o timer_debug.o
	$(CXX) test_matching_debug.o flow_debug.o dynamictree_debug.o timer_debug.o -o test_matching_debug
	@echo "debug build links"

.PHONY: test bench debug

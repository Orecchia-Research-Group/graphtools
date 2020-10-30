% MATLAB FUNCTION: TimeTest
%
% PURPOSE: Calculate the difference and percentage of runtime comparing the runtimes of Pairing when running with and without matching.
%
% INPUTS:
%    (string) graphFileName - name of file containing graph
%    (string) ptnFileName - name of file contating starting partition/bisection
%    (double) weirdrat_num - starting weirdrat numerator
%    (int64) weirdrat_den - starting weirdrat denominator
%
% OUTPUTS: 
%    (uint64) tNoMatch - runtime of Pairing without matching
%    (uint64) tMatch - runtime of Pairing with matching
%    (uint64) tDiff - difference between runtimes of Pairing with and without matching
%    (double) tPercent - percent that runtime of Pairing with matching is slower than without matching
%
% VARIABLES:
%    vector (int64) partitions - column vector that is the starting partition
%    sparse matrix (double) matching - matching routed
%    sparse matrix (double) G - instance graph
%    (int64) flow - flow output by Pairing
%    vector (int64) cut - mincut found by Pairing
%    tStart (uint64) - start time of Pairing
%

function [] = TimeTest(graphFileName, ptnFileName, weirdrat_num, weirdrat_den)

%INTITIALIZATION OF VARIABLES
partitions = int64(readPtn(ptnFilename));
G = int64(loadeg2graph(graphFileName));
cap_add = int64(weirdrat_num);
cap_orig = weirdrat_den;

%RUNNING CUTFIND TIMINGS WITHOUT MATCHING
tStart = (tic);
[flow, cut] = Pairing(G, partitions, cap_add, cap_orig);
tEnd = toc(tStart);
tNoMatch = int64(tEnd);

%RUNNING CUTFIND TIMINGS WITH MATCHING
tStart = tic;
[flow, cut, matching] = Pairing(G, partitions, cap_add, cap_orig);
tEnd = toc(tStart);
tMatch = int64(tEnd);

%PRINT CUTFIND TIMINGS
tPercent = 100*(tMatch-tNoMatch)/tMatch;
fprintf("Runtime of Pairing without matching %u\n",tMatch);
fprintf("Runtime of Pairing with matching %u\n", tNoMatch);
fprintf("Percent that runtime of matching is slower than without matching %d\n", tPercent);

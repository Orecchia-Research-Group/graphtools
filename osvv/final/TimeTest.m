% MATLAB FUNCTION: TimeTest
%
% PURPOSE: Calculate the difference and percentage of runtime comparing the runtimes of Pairing when running with and without matching.
%
% INPUTS:
%    (string) graphFileName - name of file containing graph
%    (string) ptnFileName - name of file contating starting partition/bisection
%    (double) weirdrat_num - starting weirdrat numerator
%    (int64) weirdrat_den - starting weirdrat denominator
%    (double) runNumber - number of times both Pairings are run
%
%
% OUTPUTS: 
%    (double) tNoMatch - runtime of Pairing without matching
%    (double) tMatch - runtime of Pairing with matching
%    (double) tDiff - difference between runtimes of Pairing with and without matching
%    (double) tPercent - percent that runtime of Pairing with matching is slower than without matching
%    (double) avgNM - average runtime without matching
%    (double) stdNM - standard deviation without matching
%    (double) avgM - average runtime with matching
%    (double) stdM - standard deviation with matching
%
% VARIABLES:
%    vector (int64) partitions - column vector that is the starting partition
%    sparse matrix (double) matching - matching routed
%    sparse matrix (double) G - instance graph
%    (int64) flow - flow output by Pairing
%    vector (int64) cut - mincut found by Pairing
%    (uint64) tStart - start time of Pairing
%    matrix (double) runTimeArr - matrix of runtimes
%

function [avgInit, avgS1, avgS2, avgMatch] = TimeTest(graphFileName, ptnFileName, weirdrat_num, weirdrat_den, runNumber, matching_algorithm)
if(~exist("runNumber", "var"))
	runNumber = 1;
end

%INTITIALIZATION OF VARIABLES
partitions = readPtn(ptnFileName);
partition = int64(partitions{1});
[G, n, m] = loadeg2graph(graphFileName);
[cap_add, cap_orig, minweirdrat, ex_num, ex_den, ex, cut, matching, matchrat, iterflownumber] =  RunFlow(G, partition, int64(10), int64(1), 10, int64(10), 0, matching_algorithm);
[flow, cut] = Pairing(G, partition, cap_add, cap_orig);
runTimeArr = zeros(runNumber, 4);

for i=1:runNumber

    %RUNNING CUTFIND TIMINGS WITH MATCHING
        
        [flow, cut, matching, t_init, t_S1, t_S2, t_match] = Pairing(G, partition, cap_add, cap_orig, matching_algorithm);
    
    runTimeArr(i, 1) = t_init;
    runTimeArr(i, 2) = t_S1;
    runTimeArr(i, 3) = t_S2;
    runTimeArr(i, 4) = t_match;
end

avgInit = mean(runTimeArr(1:end, 1));
avgS1 = mean(runTimeArr(1:end, 2));
avgS2 = mean(runTimeArr(1:end, 3));
avgMatch = mean(runTimeArr(1:end, 4));
end

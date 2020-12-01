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

function [avgNM, avgM, stdNM, stdM] = TimeTest(graphFileName, ptnFileName, weirdrat_num, weirdrat_den, runNumber)

if(~exist("runNumber", "var"))
	runNumber = 1;
end

%INTITIALIZATION OF VARIABLES
partitions = readPtn(ptnFileName);
partition = int64(partitions{1});
[G, n, m] = loadeg2graph(graphFileName);
[cap_add, cap_orig, minweirdrat, ex_num, ex_den, ex, cut, matching, matchrat, iterflownumber] =  RunFlow(G, partition, int64(10), int64(1), 10, int64(10), 0);
[flow, cut] = Pairing(G, partition, cap_add, cap_orig);
runTimeArr = zeros(runNumber, 2);

for i=1:runNumber
    %RUNNING CUTFIND TIMINGS WITHOUT MATCHING
    tStart = (tic);
    [flow, cut] = Pairing(G, partition, cap_add, cap_orig);
    tNoMatch = toc(tStart);

    %RUNNING CUTFIND TIMINGS WITH MATCHING
    tStart = tic;
    [flow, cut, matching] = Pairing(G, partition, cap_add, cap_orig);
    fprintf('flow value is equal to %f\n', flow);
    tMatch = toc(tStart);
    
    runTimeArr(i, 1) = tNoMatch;
    runTimeArr(i, 2) = tMatch;
end

avgNM = mean(runTimeArr(1:end, 1));
stdNM = std(runTimeArr(1:end, 1));
avgM = mean(runTimeArr(1:end, 2));
stdM = std(runTimeArr(1:end, 2));

%PRINT CUTFIND TIMINGS
%tPercent = 100*(tMatch-tNoMatch)/tMatch;
%fprintf('Runtime of Pairing without matching %.3fs\n',tNoMatch);
%fprintf('Runtime of Pairing with matching %.3fs\n', tMatch);
%fprintf('Percent that runtime of matching is slower than without matching %.2f%%\n', tPercent);
end

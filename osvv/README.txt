SETUP FOR RUNNING EXPERIMENTS

CODE ORGANIZATION
THe code is organized in 4 directories: final, logs, results and scripts. There is also a directory data containing examples.
The FINAL directory contains the matlab and C source code, together with the makefile to compile the code on your system.
The LOGS direcotry is used to store the log files of the runs of the algorithm. Similarly, the RESULTS directory is used to store the result files.
For each graph and setting of the algorithm parameters there exists a log file and a result file to which the algorithms appends its new log and results.
The SCRIPTS directory contains bash scripts that facilitate the usage of the cutfind program.


CHECKLST
INSTALLATION:
- Replace first line of bash scripts in SCRIPTS directory with address of BASH shell on your system.
- Set environement variable MATLAB to your MATLAB-r2007a directory. Also, the command "matlab" should run matlab-r2007a from the shell.
- Go to FINAL and run make.
- Test the correct functioning: Go to SCRIPTS and type "fullsuite ../data/small-graqphs football 0 5". This will run the full suite of eperiments on the small graph football.eg2. 
It should take a few minutes and afterwards you should see a bunch of files in the LOGS and RESULTS directory containing the results. You can use statsscript and bestof to compute statistics
by running them on the result files in the RESULTS directory.


PARAMETERS:
The parameters of the algorithm are:
- max number of iterations
- stopping condition
- learning rate
- initialization coefficent (relative weight of isntance graph compared to matchings added)
- flow precision
- rate spec
- lower bound spec
- certificate spec
For our purpose only the top 4 matter, as we will have flow precision always set to 512, rate_spec always set to 'n' (usual learning rate) and lower bound set to 'n', as we don't want to compute lower bounds in these runs. Certificate spec is set to 1 as we are not interested in the algorithm returning the certificate of expansion.

RESULT FILE:
A result file's filename indicates which graph and parameter choice the information refers to. For example,
geo10kplus.000.100.10.10.1.512.'n'.'n' A result file is composed of many lines, each for a run.
a typical line looks like:
r:      83      0.204579        697     3407    1.383395        0.678818        1       0.000000        6       0.575590        0.796544        0.000685
r: indictes that the line refers to the result of a run.
83 is the random seed used.
0.204579 is the score obtained
697 is the numerator of the score
3407 is the denominator
1.383395 is the time taken (without init)
0.6 is the init time
1 is the number of iterations actually run before stopping condition kicked in
0.0000 is the lower bound value obtained. This is non-zero when the lowerbound spec is set to 'y' or 'ylast'.
6 is the number of maxflow computations run
0.575590 is the time used by spectral computations
0.796544 is the time used by flow computations
0.00685 is the time spent computing the lowerbound. This may be non-zero even if a lower bound is not required as it includes the evaluation of an if clause.


FULLSUITE
Fullsuite has been modified to take 6 inputs.
fullsuite graph_location graph_name(no eg2 extension) starting_seed #trials starting_seed_special #trials_special
The special parameters refer to special runs of the algorithms: these are one shot runs, i.e runs with no feedback, where the exponential walk is computed once
and a sodaimprovement computation is then performed. There results will be a good test of how much the feedback helps and of whether the approximate eigenvector computation
given by the exponential can imrpove over specsoda in terms of time or score for some graphs, especially GMiller.


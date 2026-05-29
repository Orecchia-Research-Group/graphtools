% MATLAB FUNCTION: cutfind
%
% PURPOSE: runs algorithm OSVV on graph
%
% INPUTS:
%    (char) fileToRead - a eg2 undirected graph file to read - must be a valid graph
%    (char | 1 | 2) outputFile - output file for run results - stopping condition is appended
%    (char) suffix - Output filename suffix
%    (int16) t - maximum number of iterations of algorithms - positive - max 2^16 -1 - rounded if not integral
%    (double) stop - number of iterations after which, if no improv. in weirdratio, program exits - rounded below if not integral - can be array now
%    (double) eta - learning parameter - must be positive
%    (double) init - specification of weight of G in initialization
%    (int32) seed - seed for random number generator
%    (int64) p - precision used to compute expander flows - positive -  rounded if not integral
%    (int64) pwr_k - How many cut vectors to use in each iteration
%    (char) rate - specification of the learning rate to be used
%                                     'd' - equals eta sqrt(8log(n)/t)
%                                     'infty' - uses the second smallest eigenvalue of the Laplacian;
%                                     'n'  - equals eta;
%    (char) lwbd - 'ylast' if final lower bound desired.
%                  'y' if lower bound at 'stop' is desired.
%                  'yall' if lower bound at all iterations is desired.
%                  'n' otherwise
%    (char) matchingAlgorithm - algorithm to use for flow decomposition
%                                     'dinic' - start from source walk to sink; start again
%                                     'dynamic' - Use dynamic trees
%    (double) certificateSpec - 1 if certificate is required; 0 otherwise
%    (double) ufactor - fraction of total volume in smaller 
%    (int64) lambda_num - numerator for lambda controlling how much flow can pass through a node
%    (int64) lambda_den - denominator for lambda controoling how much flow can pass through a node

% OUTPUTS:
%    expansionFound - best  expansion score found
%    edgesCut - number of edges in best cut found
%    communities - best communities found
%    H - certificate of expansion, a graph on same vertex set of G
%    endtime - total time taken
%    inittime - time taken by initialization
%    spectime - time taken by spectral computation
%    flowtime - time taken by flow computation
%    iterations - total number of iterations run
%    flownumber - total number of maxflow computations run
%    scores - score at each iteration
%    iterscores - tx3 with iteration number, score and lower bound at stop intervals

% ISSUES:
% - should have max number of vertices or edges?)
% - nmin label is assumed to be 0 or 1? should be 0 outside the program. 1 in matlab
% - does it make sense to use weirdrat as bound guiding the search?

function [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, iterscores, lower] = ...
    cutfind(fileToRead, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ARGUMENTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    fileToRead (:, :) {mustBeFileOrGraph}
    options.outputfile (1, :) {mustBeFileOrID(options.outputfile, 0, 1)} = 1
    options.suffix (1, :) char = ''
    options.t (1, 1) int16 {mustBePositive} = 100
    options.stop (1, :) int64 {mustBeNumeric, mustBeGreaterThanOrEqual(options.stop, 1)} = 10
    options.eta (1, 1) double {mustBePositive} = 0.5
    options.init (1, 1) double {mustBeNonnegative} = 1
    options.seed (1, 1) {mustBeNumeric} = 0
    options.p (1, 1) int64 {mustBeGreaterThanOrEqual(options.p, 1)} = 1000
    options.pwr_k (1, 1) int64 {mustBePositive} = 1
    options.rate (1, :) char {mustBeMember(options.rate, {'d', 'n', 'infty', 'KL'})} = 'n'
    options.lwbd (1, :) char {mustBeMember(options.lwbd, {'ylast', 'y', 'yall', 'n'})} = 'y'
    options.matchingAlgorithm (1, :) char {mustBeMember(options.matchingAlgorithm, {'dinic', 'dynamic'})} = 'dinic'
    options.certificateSpec (1, 1) {mustBeNumericOrLogical, mustBeInRange(options.certificateSpec, 0, 1)} = 1
    options.ufactor (1, 1) double {mustBeLessThanOrEqual(options.ufactor, 0.5)} = 0
    options.lambda_num (1, 1) int64 = 1
    options.lambda_den (1, 1) int64 {mustBePositive} = 1
    options.verbose (1, 1) int64 {mustBeNonnegative} = 0
end

outputfile = options.outputfile;
suffix = options.suffix;
t = int16(options.t);
stop = options.stop;
eta = options.eta;
init = options.init;
seed = options.seed;
p = int64(options.p);
pwr_k = options.pwr_k;
rate = options.rate;
lwbd = options.lwbd;
matchingAlgorithm = options.matchingAlgorithm;
certificateSpec = options.certificateSpec;
ufactor = options.ufactor;
lambda_num = options.lambda_num;
lambda_den = options.lambda_den;
verbose = options.verbose;

size_stop = size(stop,2);

%%%%%%%%%%%%%%%%%%%%%% READ GRAPH & INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%

% TIMER ON
tic;

% RANDOM NUMBER GENERATOR INITIALIZATION
rand('twister', seed);

% READ GRAPH G
[G, weight] = loadGraph(fileToRead);
n = size(G, 1);

% Check if G is directed or not
isDirected = (nnz(G - G') ~= 0);

% CONVERT m FROM NUMBER OF ARCS TO NUMBER OF EDGES
m = nnz(G)/2;
degree = int64(full(sum(G)));
sparse_deg = diag(sum(G));
vol = sum(weight);
factor = sparse(1:n, 1:n, double(weight).^(-1/2));

%  INTIAL CERTIFICATE
H = init*G;
D = diag(sum(H));

% INITIALIZE EXPANSION AND WEIRDEST RATIO TRACKER VARIABLES - DOES NOT WORK FOR BALANCED
% MINEXP
[~, bestcut] = max(sum(G)); % find maximum degree and maximum degree vertex
minexp_num = max(full(sum(G)));
bestcut = int64(bestcut);

minexp_den = int64(1);
minexp = minexp_num;

%MINWEIRDRATfactor = diag(double(weight).^(-1/2));

if(minexp_num > p)
    fprintf(2,'Max degree is higher than precision. Search will start at weirdrat = p.\n');
    
    minweirdrat = double(p);
    minweirdrat_num = double(p);
    minweirdrat_den = int64(1);
    
else
    
    minweirdrat = minexp;
    minweirdrat_num = minexp_num;
    minweirdrat_den = int64(minexp_den);
    
end

%STOPPING CONDITION
stop = sort(stop);
stop_cnt = 1;
notimproved = 0;
if(ischar(outputfile))
    for k=1:size_stop
        infix = sprintf('.%d.', stop(k));
        output(k) = fopen(strcat(outputfile, infix, suffix),  'at');
    end
else
    for k=1:size_stop
        output(k) = outputfile;
    end
end

% INITIAL ENTROPY
entr = log2(n);
%ppool = gcp();
%pwr_k = min(pwr_k, ppool.NumWorkers);

% LOWERBOUND
congestion = init;

%COUNTER
flownumber = 0;

% CUMULATIVE TIMERS
spectime = 0;
flowtime = 0;
lowertime = 0;

% CERTIFICATE SPECIFICATION
nomatching = 0;
u = ones(n, pwr_k);

% PARALLEL PLACEHOLDERS
matchrat = zeros(pwr_k, 1);
matching = cell(1, pwr_k);
degree_distortion = zeros(pwr_k, 1);
iterflownumber = zeros(pwr_k, 1);
cut = cell(1, pwr_k);
reciprocalCut = cell(1, pwr_k);
ex = zeros(pwr_k, 1);
ex_num = zeros(pwr_k, 1);
ex_den = zeros(pwr_k, 1);
weirdrat = zeros(pwr_k, 1);
weirdrat_num = zeros(pwr_k, 1);
weirdrat_den = zeros(pwr_k, 1);

inittime = toc;

%%%%%%%%%%%%%%%%%%%%%%% POST INITIALIZATION SUMMARY %%%%%%%%%%%%%%%%%%%%%%%
if(verbose > 0)
    fprintf(2, '\nInitialization complete. Time required: %f\n', inittime);
    fprintf(2, '\nRunning on ...\n');
    fprintf(2, 'Number of vertices: %d. Number of edges: %d. Graph volume: %d\n', n, m, vol);
    fprintf(2, 'Number of iterations: %d.\n', t);
    
    fprintf(2, 'Stopping condition:');
    for k=1:size_stop
        fprintf(2,' %d', stop(k));
    end
    fprintf(2,'\n');
    
    fprintf(2, 'Learning rate: %f.\n', eta);
    fprintf(2, 'Initialization: %f.\n', init);
    fprintf(2, 'Random generator seed: %f.\n', seed);
    fprintf(2, 'Flow precision: %ld.\n', p);
    fprintf(2, 'Run rate: %s.\n', rate);
    fprintf(2, 'Lower bound: %s.\n', lwbd);
    fprintf(2, 'Vector number: %d.\n', pwr_k);
    fprintf(2, 'Matching algorithm: %s\n', matchingAlgorithm);
    fprintf(2, 'Lambda: %d / %d = %.2f.\n', lambda_num, lambda_den, double(lambda_num) / double(lambda_den));
end
%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM MAIN LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
for i=1:double(t)
    %% Cut Player Algorithm 3, pp 30 of Chen et al
    tSpectral = tic;
    [A, B] = cutPlayer(H, i, weight, factor, sparse_deg, init=init, rate=rate, eta=eta, pwr_k=pwr_k, embedding_dim=1, directed=isDirected, t=3, b=1/10);
    spectime = spectime + toc(tSpectral);

    %% Parallel vector cut/matching
    parfor step=1:pwr_k
        nomatching = 0;
        % IF CERTIFICATESPEC = 1 DO NOT NEED TO COMPUTE MATCHING IN LAST ITERATION - USED ESPECIALLY in NO FEEDBACK RUNS
        if(strcmp(lwbd,'n') && certificateSpec == 1 && i == t)
            nomatching = 1;
        end

        tFlow = tic;
        % CALL SODA_IMPROV AND ROUTING PROCEDURE IN RUNFLOW
        [weirdrat_num(step), weirdrat_den(step), weirdrat(step), ex_num(step), ex_den(step), ex(step), ...
            cut{step}, reciprocalCut{step}, matching{step}, matchrat(step), iterflownumber(step)] = ...
            RunFlow(G, A, B, weight, minweirdrat_num, minweirdrat_den, minweirdrat, p=p, ...
            nomatching_flag=nomatching, matching_algorithm=matchingAlgorithm, ufactor=ufactor, ...
            lambda_num=lambda_num, lambda_den=lambda_den);
        flowtime = flowtime + toc(tFlow);
        % fprintf(1, "%d %d\n", nnz(matching), size(matching, 2));
        % UPDATE CERTIFICATE
    
        % fprintf(2, 'Min = %d Max = %d\n', full(min(sum(matching))), full(max(sum(matching))));
        degree_distortion(step) = full(max(double(sum(matching{step})) ./ double(weight)));
        if(verbose > 1)
            fprintf(2, 'Metric = %f\n', full(max(double(sum(matching{step})) ./ double(weight))));
            fprintf(2, 'Nonzero element of matching: %d. Nonzero elements of sum %d |matching|_inf = %f\n', nnz(matching{step}), nnz(H), norm(factor * matching{step} ./ degree_distortion(step) * factor, inf));
        end
    end
    
    %% Update from parallel
    % UPDATE LOWER BOUND
    for step=1:pwr_k
        u_factor(step) = 1 / double(pwr_k);
        congestion = congestion + 1 / matchrat(step) * u_factor(step);
        % fprintf(2, 'Volume of matching %d\n', sum(matching, 'all'));
        H = H + double(matching{step}) *  u_factor(step) ./ degree_distortion(step);
        D = D + double(diag(sum(matching{step}))) *  u_factor(step) ./ degree_distortion(step);
        % UPDATE COUNTER
        flownumber = flownumber + iterflownumber(step);
    end
    % CHECK IF CUT FOUND BEATS BEST CUT
    if(~isempty(cut)) % if some cut has been found
        [minexp, minindex] = min([ex; minexp]);
        if minindex < pwr_k + 1
            bestcut=cut{minindex};
            reciprocalBestcut = reciprocalCut{minindex};
            minexp = ex(minindex);
            minexp_num = ex_num(minindex);
            minexp_den = ex_den(minindex);
            notimproved = 0;
            minweirdrat = weirdrat(minindex);
            minweirdrat_num = weirdrat_num(minindex);
            minweirdrat_den = weirdrat_den(minindex);
        else
            notimproved = notimproved + 1;
        end
        if(~strcmp(lwbd, 'n'))
            certificate = factor * ((init + i - 1) .* sparse_deg - H) * factor;
            certificatecongestion = congestion;
        end
    else
        notimproved = notimproved + 1;
    end
    % PRINT CURRENT RESULT
    if(verbose > 0)
        fprintf(2, 'Wrat: %f. Iter %d. Exp: %d / %d = %f\n', minweirdrat, i, minexp_num, minexp_den, minexp);
    end
        
    % CHECK STOPPING CONDITION
    if(notimproved >= stop(stop_cnt) || i == t || strcmp(lwbd, 'yall'))
        endtime = toc;
        if(verbose > 0)
            fprintf(2,'\nBest cut found: %d / %ld. Expansion: %f.\n', minexp_num, minexp_den, minexp);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%% LOWER BOUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tLower = tic;
        if(strcmp(lwbd, 'y'))
            certificate = factor * ((init + i - 1) .* sparse_deg - H) * factor;
            certificatecongestion = congestion;
        end
        
        if(~strcmp(lwbd,'n'))
            opts.k = 1;
            opts.tol = 0.01;
            opts.sigma = 'se';
            [temp, eigen] = eigs(@(x) ((D - H) * x + sum(sparse_deg * x) * sparse_deg * ones(size(x))), n, sparse_deg, 1, 'SA', opts); % irbleigs(certificate, opts);
            lower = 0.5*eigen(1)/certificatecongestion;
        else
            lower = 0;
        end
        iterscores(stop_cnt, :) = [i, stop(stop_cnt), minexp, lower, endtime];
        lowertime = lowertime + toc(tLower);
        
        % PRINT RUN RESULTS TO OUTPUT FILE
        % fprintf(output(stop_cnt), 'r:\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'seed', 'minexp', 'minexp_num', 'minexp_den', 'endtime', 'inittime', 'i', 'lower', 'flownumber', 'spectime', 'flowtime', 'lowertime');
        fprintf(output(stop_cnt), 'r:\t%d\t%f\t%d\t%d\t%f\t%f\t%d\t%.8f\t%d\t%f\t%f\t%f\t%d\n', seed, minexp, minexp_num, minexp_den, endtime, inittime, i, lower, flownumber, spectime, flowtime, lowertime, nnz(H));
        
        if(~strcmp(lwbd, 'n') && verbose > 0)
            fprintf(2,'Lower bound: %f.\n', lower);
        end
        if(verbose > 0)
            fprintf(2,'Algorithm has completed. Time required: %f\n', endtime);
        end
        
        if(notimproved >= stop(stop_cnt))
            stop_cnt = stop_cnt + 1;
        end
        if(stop_cnt > size_stop)
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% TERMINATION & OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% SET OUTPUT VARIABLES
% EXPANSION
expansionFound = minexp;
edgesCut = minexp_num;
if(size(bestcut,1) < double(n)/2)
    L = bestcut;
    R = reciprocalBestcut;
else
    L = reciprocalBestcut;
    R = bestcut;
end

% ITERATIONS
iterations = i;



end

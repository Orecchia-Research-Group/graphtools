% MATLAB FUNCTION: cutfind
%
% PURPOSE: runs algorithm OSVV on graph
%
% INPUTS:
%    (char) FileToRead - a eg2 undirected graph file to read - must be a valid graph
%    (int16) t - maximum number of iterations of algorithms - positive - max 2^16 -1 - rounded if not integral
%    (double) stop - number of iterations after which, if no improv. in weirdratio, program exits - rounded below if not integral - can be array now
%    (double) eta - learning parameter - must be positive
%    (double) init - specification of weight of G in initialization
%    (int64) p - precision used to compute expander flows - positive -  rounded if not integral
%    (int64) pwr_k - How many cut vectors to use in each iteration
%    (int32) seed - seed for random number generator
%    (char) rate - specification of the learning rate to be used
%                                     'd' - equals eta sqrt(8log(n)/t)
%                                     'infty' - uses the second smallest eigenvalue of the Laplacian;
%                                     'n'  - equals eta;
%    (int32) lambda - flow that can pass through a node. If missing considered to be infinity


%    (char) lwbd - 'y' if final lower bound desired. 'n' otherwise
%    (char) outputfile - output file for run results - stopping condition is appended
%    (double) certificate - 1 if certificate is required; 0 otherwise
%    (char) suffix - other string to be appended to output file

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
% - reorder parameters

function [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, iterscores, lower] = ...
    cutfind(FileToRead, outputfile, suffix, t, stop,  eta, init, seed, p, pwr_k, rate, lwbd, certificatespec, ufactor, varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%  ERROR CHECKING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_string = 'Error in parameter %s. See README file for usage.\n';

if(~ischar(FileToRead) && (~isnumeric(FileToRead) || (size(FileToRead, 1) ~= size(FileToRead, 2))))
    error(error_string, 'eg2 input file');
end

if(~isnumeric(t) || t < 1)
    error(error_string, 'number of iterations');
end
t = int16(t);

if(~isnumeric(eta) || eta < 0)
    error(error_string, 'learning rate');
end

if(~isnumeric(p) || p < 1)
    error(error_string, 'flow precision');
end
p = int64(p); %rounded to integer - should check?

if(~isnumeric(pwr_k) || pwr_k < 1)
    error(error_string, 'vector number');
end
pwr_k = double(pwr_k);

if(~isnumeric(seed))
    error(error_string, 'seed');
end

if(~ischar(rate)|| ~(strcmp(rate, 'd') || strcmp(rate, 'n') || strcmp(rate,'infty') || strcmp(rate,'KL')))
    error(error_string, 'rate spec');
end

if(~isnumeric(init) || init < 0)
    error(error_string, 'init weight');
end

if(~isnumeric(stop))
    error(error_string, 'stopping condition');
end
size_stop = size(stop,2);
if(size_stop <1)
    error(error_string, 'stopping condition');
end

for k=1:size_stop
    if(stop(k) < 1)
        error(error_string, 'stopping condition');
    end
end

if(~ischar(lwbd) || ~(strcmp(lwbd, 'y') || (strcmp(lwbd, 'n') || strcmp(lwbd, 'ylast'))))
    error(error_string, 'lowerbound computation');
end

if(~ischar(outputfile) && outputfile ~= 1 && outputfile ~= 2)
    error(error_string, 'output file name');
end

if(certificatespec ~= 1 && certificatespec ~= 0)
    error(error_string, 'certificate specification');
end


if(~ischar(suffix))
    error(error_string, 'suffix file name');
end

% Mixed cut or edge cut?
if (size(varargin, 2) > 0)
    if length(varargin) < 2
        error(error_string, 'lambda');
    end
    lamda_num = int64(varargin{1});
    lamda_den = int64(varargin{2});
    if lamda_num > lamda_den
        error('Lambda needs to be less than or equal to 1');
    end
else
    lamda_num = int64(-1);
    lamda_den = int64(1);
end


%%%%%%%%%%%%%%%%%%%%%% READ GRAPH & INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%

% TIMER ON
tic;

% RANDOM NUMBER GENERATOR INITIALIZATION
rand('twister', seed);

% READ GRAPH G
if(ischar(FileToRead))
    [G, weight] = loadMetisGraph(FileToRead);
    n = size(G, 1);
    weight = int64(weight);
else
    G = FileToRead;
    n = size(G, 1);
    degree = int64(full(sum(G)));
    weight = ones(1, n, 'int64');
    weight(:) = degree;
end
% ERROR CHECKING: G MUST BE UNDIRECTED
% if(nnz(G - G') ~= 0)
%    error('The eg2 graph is not undirected.\n');
% end

% CONVERT m FROM NUMBER OF ARCS TO NUMBER OF EDGES
m = nnz(G)/2;
degree = int64(full(sum(G)));
sparse_deg = diag(sum(G));
vol = sum(weight);
factor = diag(sum(G).^(-1/2));

% Check that the graph is strongly connected
[S, C] = graphconncomp(G);
if (S > 1)
    %[GC, GR] = groupcounts(C);
    %[out,idx] = sort(GC);
    error('Graph is not strongly connected. S = %d', S);
end

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

%MINWEIRDRAT
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

%%%%%%%%%%%%%%%%%%%%%%% POST INITIALIZATION SUMMARY %%%%%%%%%%%%%%%%%%%%%%%
inittime = toc;
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
fprintf(2, 'Lambda: %d / %d = %.2f.\n', lamda_num, lamda_den, double(lamda_num) / double(lamda_den));

%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM MAIN LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
for i=1:double(t)
    
    tSpectral = tic;
    % RANDOM BISECTION INITIALIZATION;
    v = round(rand(n, pwr_k));
    v = v - mean(v);
    
    
    % LEARNING RATE INITIALIZATION
    if(strcmp(rate,'d'))
        current_eta = eta*sqrt(8*log(n)/i);
    elseif (strcmp(rate, 'KL'))
        current_eta = eta*sqrt(8*entr/i);
    else    
        current_eta = eta;
    end
    
    % SPECTRAL PARTITIONING
    %% SECOND EIGENVALUE
    M = factor * ((init + i - 1) .* sparse_deg - H) * factor;
    if(strcmp(rate,'infty'))
        % opts.k = 2;
        opts.tol = 1e-14;
        % opts.sigma = 'SM';
        [u, ~] = eigs(@(x) (((init + i - 1) .* sparse_deg - H) * x + sum(sparse_deg * x) * sparse_deg * ones(size(x))), n, sparse_deg, pwr_k, 'SA', opts);
        u = factor * u;
        u_factor(1:pwr_k) = 1 / pwr_k;
    end
    %% Parallel vector cut0matching
    parfor step=1:pwr_k
        nomatching = 0;
        if(~strcmp(rate, 'infty'))
            %%%  RANDOM WALK STEP %%% RESCALE V for better tolerance
            u(:, step) = factor * expv((-1)*current_eta, M, v(:, step), 1e-3);
            u_factor(step) = 1 / pwr_k;
        end
        spectime = spectime + toc(tSpectral);


        % SORT VERTICES BY PROB. CHARGE AND DETERMINE BISECTION
        [~, index] = sort(u(:, step));
        index = int64(index);
        j = floor(n / 2) + 1;
        bisec=int64(index(1:floor(n/2)));
        bisec_vol = sum(weight(bisec));

        while bisec_vol < vol / 2.0
            bisec(end + 1) = int64(index(j));
            bisec_vol = bisec_vol + weight(index(j));
            j = j + 1;
        end
        while bisec_vol > vol / 2.0
            bisec_vol = bisec_vol - weight(bisec(end));
            bisec = bisec(1:end-1);
        end
        fprintf(2, 'Bisec volume: %ld. Bisec size: %ld. Vol frac: %f.\n', full(bisec_vol), length(bisec), full(double(bisec_vol) / vol));

        % IF CERTIFICATESPEC = 1 DO NOT NEED TO COMPUTE MATCHING IN LAST ITERATION - USED ESPECIALLY in NO FEEDBACK RUNS
        if(strcmp(lwbd,'n') && certificatespec == 1 && i == t)
            nomatching = 1;
        end

        tFlow = tic;
        % CALL SODA_IMPROV AND ROUTING PROCEDURE IN RUNFLOW
        if (lamda_num > 0)
            [weirdrat_num(step), weirdrat_den(step), weirdrat(step), ex_num(step), ex_den(step), ex(step), cut{step}, reciprocalCut{step}, matching{step}, matchrat(step), iterflownumber(step)] =  ...
                RunFlow(G, bisec, weight, minweirdrat_num, minweirdrat_den, minweirdrat, p, nomatching, ufactor, lamda_num, lamda_den);
        else
            [weirdrat_num(step), weirdrat_den(step), weirdrat(step), ex_num(step), ex_den(step), ex(step), cut{step}, reciprocalCut{step}, matching{step}, matchrat(step), iterflownumber(step)] =  ...
                RunFlow(G, bisec, weight, minweirdrat_num, minweirdrat_den, minweirdrat, p, nomatching, ufactor);
        end
        flowtime = flowtime + toc(tFlow);
        % fprintf(1, "%d %d\n", nnz(matching), size(matching, 2));
        % UPDATE CERTIFICATE
    
        % fprintf(2, 'Min = %d Max = %d\n', full(min(sum(matching))), full(max(sum(matching))));
        fprintf(2, 'Metric = %f\n', full(max(double(sum(matching{step})) ./ double(weight))));
        degree_distortion(step) = full(max(double(sum(matching{step})) ./ double(weight)));
        fprintf(2, 'Nonzero element of matching: %d. Nonzero elements of sum %d |matching|_inf = %f\n', nnz(matching{step}), nnz(H), norm(factor * matching{step} ./ degree_distortion(step) * factor, inf));
    end
    
    %% Update from parallel
    % UPDATE LOWER BOUND
    for step=1:pwr_k
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
        if(strcmp(lwbd, 'ylast'))
            certificate = factor * ((init + i - 1) .* sparse_deg - H) * factor;
            certificatecongestion = congestion;
        end
    else
        notimproved = notimproved + 1;
    end
    % PRINT CURRENT RESULT
    fprintf(2, 'Wrat: %f. Iter %d. Exp: %d / %d = %f. eta: %f\n', minweirdrat, i, minexp_num, minexp_den, minexp, current_eta);
        
    % CHECK STOPPING CONDITION
    if(notimproved >= stop(stop_cnt) || i == t)
        endtime = toc;
        fprintf(2,'\nBest cut found: %d / %ld. Expansion: %f.\n', minexp_num, minexp_den, minexp);
        
        
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
        lowertime = toc(tLower);
        
        % PRINT RUN RESULTS TO OUTPUT FILE
        % fprintf(output(stop_cnt), 'r:\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'seed', 'minexp', 'minexp_num', 'minexp_den', 'endtime', 'inittime', 'i', 'lower', 'flownumber', 'spectime', 'flowtime', 'lowertime');
        fprintf(output(stop_cnt), 'r:\t%d\t%f\t%d\t%d\t%f\t%f\t%d\t%f\t%d\t%f\t%f\t%f\n', seed, minexp, minexp_num, minexp_den, endtime, inittime, i, lower, flownumber, spectime, flowtime, lowertime);
        
        if(~strcmp(lwbd, 'n'))
            fprintf(2,'Lower bound: %f.\n', lower);
        end
        fprintf(2,'Algorithm has completed. Time required: %f\n', endtime);
        
        
        
        stop_cnt = stop_cnt + 1;
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

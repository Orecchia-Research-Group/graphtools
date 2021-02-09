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
%    (int32) seed - seed for random number generator
%    (char) rate - specification of the learning rate to be used
%                                     'd' - equals eta sqrt(8log(n)/t)
%                                     'infty' - uses the second smallest eigenvalue of the Laplacian;            
%                                     'n'  - equals eta; 


%    (char) lwbd - 'y' if final lower bound desired. 'n' otherwise
%    (char) outputfile - output file for run results - stopping condition is appended
%    (double) certificate - 1 if certificate is required; 0 otherwise
%    (char) suffix - other string to be appended to output file

% OUTPUTS:
%    expansionFound - best  expansion score found
%    edgesCut - number of edges in best cut found
%    cutFound - list of vertices composing best cut found
%    H - certificate of expansion, a graph on same vertex set of G
%    endtime - total time taken 
%    inittime - time taken by initialization
%    spectime - time taken by spectral computation
%    flowtime - time taken by flow computation
%    iterations - total number of iterations run
%    flownumber - total number of maxflow computations run
%    lower - lower bound found

% ISSUES:
% - should have max number of vertices or edges?
% - nmin label is assumed to be 0 or 1? should be 0 outside the program. 1 in matlab
% - does it make sense to use weirdrat as bound guiding the search?
% - reorder parameters

function [expansionFound, edgesCut, cutFound,H, iterations, lower] = cutfind(FileToRead, outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, matchingAlgorithm, certificatespec) 


%%%%%%%%%%%%%%%%%%%%%%%%%%  ERROR CHECKING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_string = 'Error in parameter %s. See README file for usage.\n';

if(~ischar(FileToRead))
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

if(~isnumeric(seed))
  error(error_string, 'seed');
end
   
if(~ischar(rate)|| ~(strcmp(rate, 'd') || strcmp(rate, 'n') || strcmp(rate,'infty')))
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
   error(error_string, 'lowerbound argumentsd computation');
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

if(~ischar(matchingAlgorithm) || (~strcmp(matchingAlgorithm, 'dinic') && ~strcmp(matchingAlgorithm, 'dynamic')))
    error(error_string, 'matching algorithm');
end


%%%%%%%%%%%%%%%%%%%%%% READ GRAPH & INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%  

% TIMER ON
tic;

% RANDOM NUMBER GENERATOR INITIALIZATION
rand('twister', seed);

% READ GRAPH G
[G, n, m] = loadeg2graph(FileToRead);
 
% ERROR CHECKING: G MUST BE UNDIRECTED
if(nnz(G - G') ~= 0)
   error('The eg2 graph is not undirected.\n');
end

 % CONVERT m FROM NUMBER OF ARCS TO NUMBER OF EDGES
m = nnz(G)/2;

%  INTIAL CERTIFICATE
H = init*G;
D =diag(sum(H));

% INITIALIZE EXPANSION AND WEIRDEST RATIO TRACKER VARIABLES
% MINEXP
[minexp_num, bestcut] = max(sum(G)); % find maximum degree and maximum degree vertex
minexp_num = full(minexp_num); % make sure it's full
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

%%%%%%%%%%%%%%%%%%%%%%% POST INITIALIZATION SUMMARY %%%%%%%%%%%%%%%%%%%%%%%
inittime = toc;
fprintf(2, '\nInitialization complete. Time required: %f\n', inittime);
fprintf(2, '\nRunning on ...\n');
fprintf(2, 'Graph: %s. Number of vertices: %ld. Number of edges: %ld.\n', FileToRead, n, m);
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

%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM MAIN LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
for i=1:double(t) 
    
    tSpectral = tic;
    % RANDOM BISECTION INITIALIZATION;   
    v = round(rand(n,1));
    
    
    % LEARNING RATE INITIALIZATION
    if(strcmp(rate,'d'))      
      eta = eta*sqrt(8*log(n)/i);
    end
   
    % SPECTRAL PARTITIONING
    %% SECOOND EIGENVALUE
    if(strcmp(rate,'infty'))
       opts.k = 2;
       opts.tol = 0.001;
       opts.sigma = 'se';
       [temp, trash ] = irbleigs(D-H, opts);
       u=temp([1:n],2);
    else
%%%  RANDOM WALK STEP %%% RESCALE V for better tolerance
    u = expv((-1)*eta, (D - H), v);
%%%                   %%%
    end
    spectime = spectime + toc(tSpectral);


     % SORT VERTICES BY PROB. CHARGE AND DETERMINE BISECTION
     [ordered, index] = sort(u);
     index = int64(index);
     bisec=int64(index(1:floor(n/2)));

    % IF CERTIFICATESPEC = 1 DO NOT NEED TO COMPUTE MATCHING IN LAST ITERATION - USED ESPECIALLY in NO FEEDBACK RUNS
    if(strcmp(lwbd,'n') && certificatespec == 1 && i == t)
       nomatching = 1;
    end    
    
    tFlow = tic;
    % CALL SODA_IMPROV AND ROUTING PROCEDURE IN RUNFLOW    
    [minweirdrat_num, minweirdrat_den, minweirdrat, ex_num, ex_den, ex, cut, matching, matchrat, iterflownumber] =  RunFlow(G, bisec, minweirdrat_num, minweirdrat_den, minweirdrat, p, nomatching, matchingAlgorithm);
    flowtime = flowtime + toc(tFlow);

    % UPDATE LOWER BOUND
    congestion = congestion +1/matchrat;

    % UPDATE CERTIFICATE
    H = H + matching;
    D = D + diag(sum(matching));

   % UPDATE COUNTER 
   flownumber = flownumber + iterflownumber;

    % CHECK IF CUT FOUND BEATS BEST CUT
    if(~isempty(cut)) % if some cut has been found
       if(ex < minexp)
           bestcut=cut;
           minexp = ex;
           minexp_num = ex_num;
           minexp_den = ex_den;          
       end
       notimproved = 0;
       if(strcmp(lwbd, 'ylast'))
          certificate = D-H;
          certificatecongestion = congestion;
       end    
    else
       notimproved = notimproved + 1;
    end

    % PRINT CURRENT RESULT
    fprintf(2, 'Wrat: %f. Iter %d. Exp: %f / %d = %f\n', minweirdrat, i, minexp_num, minexp_den, minexp);


    
    % CHECK STOPPING CONDITION
   if(notimproved >= stop(stop_cnt) || i == t)
      endtime = toc;
      fprintf(2,'\nBest cut found: %f / %ld. Expansion: %f.\n', minexp_num, minexp_den, minexp);
	
   
    %%%%%%%%%%%%%%%%%%%%%%%% LOWER BOUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       tLower = tic;
       if(strcmp(lwbd, 'y'))
          certificate = D-H;
          certificatecongestion = congestion;
       end

       if(~strcmp(lwbd,'n'))
          opts.k = 2;
          opts.tol = 0.001;
          opts.sigma = 'se';
          [temp, eigen ] = irbleigs(certificate, opts);
          lower = 0.5*eigen(2,2)/certificatecongestion;    
      else
          lower = 0;
      end
      lowertime = toc(tLower);
   
   % PRINT RUN RESULTS TO OUTPUT FILE
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
   cutFound = bestcut;
else
   mask = zeros(n,1);
   mask(bestcut) = 1;
   mask = mask - ones(n,1);
   mask = (-1)*mask;
   cutFound = find(mask);
end

% ITERATIONS
iterations = i;



end

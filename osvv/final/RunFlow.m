% MATLAB FUNCTION: RunFlow
%
% PURPOSE: computes SODA improvement on a partition and returns an union of at most p matchings routed
%
% INPUTS: 
%    sparse matrix (double) G - instance graph
%    vector (int64) source_set - side conncecting to s
%    vector (int64) sink_set - side conncecting to t
%    vector (int) weight - weight of each vertex
%    (double) minweirdrat_num - starting weirdrat numerator
%    (int64) minweirdrat_den - starting weirdrat denominator
%    (double) minweirdrat - starting weirdrat
%    (int64) p - precision, specifies degree of routed matching
%    (float) lambda - flow that can pass through a node. If missing
%                     considered to be infinity
%    

%
% OUTPUTS:
%    (double) weirdrat_num - final weirdrat numerator, i.e. best weirdrat encountered in this run
%    (int64) weirdrat_den - final weirdrat denominator
%    (double) weirdrat - final weirdrat  
%    (double) ex_num - final expansion numerator
%    (int64) ex_den - final expansion denominator
%    (double) ex - final expansion, i.e. best expansion seen in this run
%    vector (int64) bestcut - best cut found in this run
%    vector (int64) reciprocalBestcut - reciprocal of the best cut found
%    sparse matrix (double) matching - matching routed   
%    (double) matchrat - weirdratio certified by matching
%    (double) flownumber - number of maxflows run
%				
% VARIABLES: 
%     (int64) n - number of vertices in G
%     (int8) found_flag - flag to terminate current run of flows
%     (int64) flow - flow output by Pairing
%     vector (int64) cut - mincut found by Pairing
%     (double) newex_num -numerator of expansion of cut found by last call of Pairing
%     (int64) newex_den - denominator of expansion of cut found by last call of Pairing
%     (double) newex -  expansion of cut found by last call of Pairing
%     (int64) cap_add - a temporary variable
%
% ISSUES: can do mincut at precision p too?
%

function   [weirdrat_num, weirdrat_den, weirdrat, ex_num, ex_den, ex, bestcut, reciprocalBestcut, matching, matchrat, flownumber] = RunFlow(G, source_set, sink_set, weight, minweirdrat_num, minweirdrat_den, minweirdrat, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%  ARGUMENTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    G (:, :) double {mustBeGraph}
    source_set (:, 1) int64
    sink_set (:, 1) int64
    weight (1, :) int64
    minweirdrat_num (1, 1) double
    minweirdrat_den (1, 1) double
    minweirdrat (1, 1) double
    options.p (1, 1) int64 = 100000
    options.nomatching_flag (1, 1) double = 1
    options.matching_algorithm (1, :) char {mustBeMember(options.matching_algorithm, {'dinic', 'dynamic'})} = 'dinic'
    options.ufactor (1, 1) double {mustBeLessThanOrEqual(options.ufactor, 0.5)} = 0
    options.lambda_num  (1, 1) int64 = -1
    options.lambda_den  (1, 1) int64 {mustBePositive} = 1
    options.verbose (1, 1) int64 {mustBeNonnegative} = 0
end

mustBeGreaterThanOrEqual(options.lambda_den, options.lambda_num)

p = options.p;
nomatching_flag = options.nomatching_flag;
matching_algorithm = options.matching_algorithm;
ufactor = options.ufactor;
lambda_num = options.lambda_num;
lambda_den = options.lambda_den;
verbose = options.verbose;

%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZATION OF OUTPUT VARIABLES
n = int64(size(G,1));
bestcut = int64([]);
reciprocalBestcut = int64([]);
matching = int64([]);
ex = Inf;
ex_num = Inf;
ex_den = int64(0);
weirdrat_num = minweirdrat_num;
weirdrat_den = int64(minweirdrat_den);
weirdrat = minweirdrat;

% PREPARE SETS
source_set = int64(sort(source_set));
source_set_size = int64(size(source_set, 1));
source_set_volume = sum(weight(source_set));

sink_set = int64(sort(sink_set));
sink_set_size = int64(size(sink_set, 1));
sink_set_volume = sum(weight(sink_set));

%%%%%%%%%%%%%%%%%%%%%%%%% SEARCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 0;

vol = sum(weight);
if (sink_set_volume > source_set_volume) && (sink_set_volume - source_set_volume) / source_set_volume < 1e-5
    side_num = int64(1);
    side_den = int64(1);
else
    [side_num, side_den] = Farey(int64(source_set_volume), int64(sink_set_volume), p);
end
weirdrat_num_lower = 0;
weirdrat_den_lower = int64(1);
weirdrat_lower = double(weirdrat_num_lower) / double(weirdrat_den_lower);

weirdrat_num_upper = weirdrat_num;
weirdrat_den_upper = int64(weirdrat_den);
weirdrat_upper = double(weirdrat_num_upper) / double(weirdrat_den_upper);

balance_violated = false;
prev_flow = -1;

while(true) % WHILE BETTER WEIRDRAT CUT EXISTS
    %%% [cap_add, cap_orig] = Farey(int64(cap_add), cap_orig, int64(10000));
    cap_add = int64(weirdrat_num);
    cap_orig = int64(weirdrat_den);
    %fprintf('weirdrat_num: %d. weirdrat_den: %d. weirdrat: %f\n', weirdrat_num, weirdrat_den, weirdrat);

    source_modifier = int64(cap_add) * int64(side_den);
    sink_modifier = int64(cap_add) * int64(side_num);
    internal_modifier = int64(cap_orig) * int64(side_den);
    original_modifier = int64(cap_orig) * int64(side_den);

    if (lambda_num > 0)
        source_modifier = int64(source_modifier * lambda_den);
        sink_modifier = int64(sink_modifier * lambda_den);
        internal_modifier = int64(internal_modifier * lambda_num);
        original_modifier = int64(original_modifier * lambda_den);
        if(verbose > 2)
            fprintf(2, 'source_modifier=%d, sink_modifier=%d, internal_modifier=%d, original_modifier=%d\n', source_modifier, sink_modifier, internal_modifier, original_modifier);
        end
        [flow, cut, reciprocalCut] = Pairing(G, source_set, sink_set, weight, matching_algorithm, source_modifier, sink_modifier, original_modifier, internal_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    else
        [flow, cut, reciprocalCut] = Pairing(G, source_set, sink_set, weight, matching_algorithm, source_modifier, sink_modifier, original_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    end
    counter = counter + 1;
    if (flow == 0)
        fprintf(2, 'You disconnected: %f\n', flow);
    end
    
    cutVolume = sum(weight(cut));
    reciprocalCutVolume = sum(weight(reciprocalCut));
    if(verbose > 2)
        fprintf(2, 'flow: %d. weirdrat_num: %d. RHS: %d. Sink side: %d. weirdrat: %f lower: %f. upper: %f\n', ...
            flow, weirdrat_num, (source_set_volume) * source_modifier, (sink_set_volume) * sink_modifier, weirdrat, weirdrat_lower, weirdrat_upper);
        fprintf(2, 'cutVolume: %d. reciprocalCutVolume: %d\n', cutVolume, reciprocalCutVolume);
    end
    
    if balance_violated
        if (min(cutVolume, reciprocalCutVolume) >= ufactor * vol)
            weirdrat_num_upper = weirdrat_num;
            weirdrat_den_upper = weirdrat_den;
            weirdrat_upper = double(weirdrat_num_upper) / double(weirdrat_den_upper);
        else
            weirdrat_num_lower = weirdrat_num;
            weirdrat_den_lower = weirdrat_den;
            weirdrat_lower = double(weirdrat_num_lower) / double(weirdrat_den_lower);
        end
    elseif(flow < double(source_set_volume) * source_modifier) % IF BETTER CUT FOUND
        %CHANGES
        [new_weirdrat_num, new_weirdrat_den, new_weirdrat] =  cutweird(G, cut, reciprocalCut, source_set, sink_set, int64(weight), int64(lambda_num), int64(lambda_den)); % COMPUTE NEW WEIRDRAT
        if new_weirdrat < weirdrat
            weirdrat_num = new_weirdrat_num;
            weirdrat_den = new_weirdrat_den;
            weirdrat = new_weirdrat;
        end

        if (min(cutVolume, reciprocalCutVolume) < ufactor * vol)
            balance_violated = true;
            weirdrat_num_lower = weirdrat_num;
            weirdrat_den_lower = weirdrat_den;
            weirdrat_lower = double(weirdrat_num_lower) / double(weirdrat_den_lower);
        else
            weirdrat_num_upper = weirdrat_num;
            weirdrat_den_upper = weirdrat_den;
            weirdrat_upper = double(weirdrat_num_upper) / double(weirdrat_den_upper);
        end

    else  % OTHERWISE STOP
        break;
    end
    if balance_violated
        rat_num = int64(weirdrat_num_lower * weird624223rat_den_upper + weirdrat_num_upper * weirdrat_den_lower);
        rat_den = 2 * weirdrat_den_lower * weirdrat_den_upper;
        [weirdrat_num, weirdrat_den] = Farey(rat_num, rat_den, p);
        weirdrat = double(weirdrat_num) / double(weirdrat_den);
    end

    [newex_num, newex_den, newex] = cutexp(G, int64(lambda_num), int64(lambda_den), int64(weight), cut, reciprocalCut);
    % CHECK IF EXPANSION HAS IMPROVED - IF IT HAS RECORD NEW CUT
    if (min(cutVolume, reciprocalCutVolume) >= ufactor * vol) && (newex < ex)
        ex_num = newex_num;
        ex_den = int64(newex_den);
        ex = newex;
        bestcut = cut;
        reciprocalBestcut = reciprocalCut;
    end
    if weirdrat_upper < (1 + 1e-3) * weirdrat_lower
        weirdrat_num = weirdrat_num_upper;
        weirdrat_den = weirdrat_den_upper;
        weirdrat = double(weirdrat_num_upper) / double(weirdrat_den_upper);
        break;
    end

    if flow == prev_flow
        break;
    end
    prev_flow = flow;
end

%% Run with upper
if(weirdrat_num == 0)
    error('RunFlow: ratio reduced to zero!');
end
if(verbose > 2)
    fprintf(2, 'Number of maxflows: %d. ', counter + 1);
end
% ONCE STOPPED, CONSTRUCT ROUTED UNION OF MATCHINGS - DO THIS AT PRECISION P
if(nomatching_flag == 0)
   %[match_num, match_den] = Farey(int64(weirdrat_num), weirdrat_den, p); % use Farey sequences to find best p-precision approximation to weirdrat
    match_num = int64(weirdrat_num);
    match_den = int64(weirdrat_den);
%    fprintf(2, 'Match_num: %f\n', match_num);
    source_modifier = int64(match_num) * int64(side_den);
    sink_modifier = int64(match_num) * int64(side_num);
    internal_modifier = int64(match_den) * int64(side_den);
    original_modifier = int64(match_den) * int64(side_den);

    if (lambda_num > 0)
        source_modifier = int64(source_modifier * lambda_den);
        sink_modifier = int64(sink_modifier * lambda_den);
        internal_modifier = int64(internal_modifier * lambda_num);
        original_modifier = int64(original_modifier * lambda_den);
        [flow, cut, reciprocalCut, matching] = Pairing(G, source_set, sink_set, weight, matching_algorithm, source_modifier, sink_modifier, original_modifier, internal_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    else
        [flow, cut, reciprocalCut, matching] = Pairing(G, source_set, sink_set, weight, matching_algorithm, source_modifier, sink_modifier, original_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    end
    matchingSum = sum(int64(full(sum(matching))));
    if 2 * flow ~= matchingSum
        fprintf(2, 'Warning: flow is %d, but matching has volume %d\n', flow, matchingSum / 2);
    end
   % MATCHING SCALING
    matching = matching / double(match_num * side_den);
    matchrat = double(match_num)/double(match_den);
    if (lambda_num > 0)
        matching = matching / double(lambda_den);
        matchrat = matchrat * double(lambda_den);
    end
    if(~issparse(matching))
        fprintf('Matching is not sparse...\n');
    end
else
    matching = sparse(1:double(n), 1:double(n), 0);
    matchrat = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT PREPARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flownumber = counter+1;

% ALL OUTPUTS ARE READY

end

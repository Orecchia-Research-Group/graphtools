% MATLAB FUNCTION: RunFlow
%
% PURPOSE: computes SODA improvement on a bisection and returns an union of at most p matchings routed
%
% INPUTS: 
%    sparse matrix (double) G - instance graph
%    vector (int64) bisec - bisection on which to run SODA improvement
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
%     (int64) size_bisec - size of the bisection bisec
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

function   [weirdrat_num, weirdrat_den, weirdrat, ex_num, ex_den, ex, bestcut, reciprocalBestcut, matching, matchrat, flownumber] = RunFlow(G, bisec, weight, minweirdrat_num, minweirdrat_den, minweirdrat, p, nomatching_flag, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FAREY_PRECISION = 10000;

if size(bisec, 2) == length(bisec)
    error('RunFlow: bisec must be a column vector\n');
end

% Mixed cut or edge cut?
if (size(varargin, 2) > 0)
    lamda = varargin{1};
else
    lamda = -1;
end
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

% PREPARE BISEC
bisec = int64(sort(bisec));
size_bisec = int64(size(bisec,1));
bisec_vol = sum(weight(bisec));

%%%%%%%%%%%%%%%%%%%%%%%%% SEARCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
found_flag = int8(1);
counter = 0;

vol = sum(weight);
if abs(2 * bisec_vol - vol) / bisec_vol < 1e-3
    side_nom = int64(1);
    side_den = int64(1);
else
    [side_nom, side_den] = Farey(int64(bisec_vol), int64(vol - bisec_vol), int64(FAREY_PRECISION));
end
%fprintf(2, 'side_nom: %d. side_den: %d.\n', side_nom, side_den);

while(found_flag) % WHILE BETTER WEIRDRAT CUT EXISTS
    %%% [cap_add, cap_orig] = Farey(int64(cap_add), cap_orig, int64(10000));
    cap_add = int64(weirdrat_num);
    cap_orig = int64(weirdrat_den);
    %fprintf('weirdrat_num: %d. weirdrat_den: %d. weirdrat: %f\n', weirdrat_num, weirdrat_den, weirdrat);
    
    source_modifier = int64(cap_add) * int64(side_den);
    sink_modifier = int64(cap_add) * int64(side_nom);
    internal_modifier = int64(cap_orig) * int64(side_den);
    original_modifier = int64(cap_orig) * int64(side_den);

    %fprintf(2, 'Source modifier %d\n', source_modifier);
    %fprintf(2, 'Sink modifier %d\n', sink_modifier);
    %fprintf(2, 'Internal modifier %d\n', internal_modifier);
    %fprintf(2, 'Original modifier %d\n', original_modifier);
    
    if (lamda > 0)
        source_modifier = int64(source_modifier / lamda);
        sink_modifier = int64(sink_modifier / lamda);
        original_modifier = int64(original_modifier / lamda);
        [flow, cut, reciprocalCut] = Pairing(G, bisec, weight, source_modifier, sink_modifier, original_modifier, internal_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    else
        [flow, cut, reciprocalCut] = Pairing(G, bisec, weight, source_modifier, sink_modifier, original_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    end
    
    % fprintf(2, 'Size of cut=%d\n', length(cut));
    
    counter = counter + 1;
    if (flow == 0)
        fprintf(2, 'You disconnected: %f\n', flow);
    end 
    
    fprintf('flow: %d. weirdrat_num: %d. RHS: %d. Sink side: %d. weirdrat: %f\n', ...
        flow, weirdrat_num, double(bisec_vol) * source_modifier, double(vol - bisec_vol) * sink_modifier, weirdrat);
    if(flow < double(bisec_vol) * source_modifier) % IF BETTER CUT FOUND
        %CHANGES
        found_flag = int8(1);
        [weirdrat_num, weirdrat_den, weirdrat] =  cutweird(G, cut, reciprocalCut, bisec, weight, lamda); % COMPUTE NEW WEIRDRAT
        % CHECK IF EXPANSION HAS IMPROVED - IF IT HAS RECORD NEW CUT
        [newex_num, newex_den, newex] = cutexp(G, lamda, weight, cut, reciprocalCut);
        if(newex < ex)
            ex_num = newex_num;
            ex_den = int64(newex_den);
            ex = newex;  
            bestcut = cut;
            reciprocalBestcut = reciprocalCut;
        end
        % PLACE HERE CODE TO PRINT EXPANSION GRAPH

    else  % OTHERWISE STOP 
        found_flag = int8(0); 
    end
end
     
	fprintf(2, 'Number of maxflows: %d. ', counter + 1);
% ONCE STOPPED, CONSTRUCT ROUTED UNION OF MATCHINGS - DO THIS AT PRECISION P
if(nomatching_flag == 0)
   [match_num, match_den] = Farey(int64(weirdrat_num), weirdrat_den, p); % use Farey sequences to find best p-precision approximation to weirdrat

%    fprintf(2, 'Match_num: %f\n', match_num);
    source_modifier = int64(match_num) * side_nom;
    sink_modifier = int64(match_num) * side_den;
    internal_modifier = int64(match_den) * side_nom;
    original_modifier = int64(match_den) * side_nom;

    if (lamda > 0)
        source_modifier = int64(source_modifier / lamda);
        sink_modifier = int64(sink_modifier / lamda);
        original_modifier = int64(original_modifier / lamda);
        [flow, cut, reciprocalCut, matching] = Pairing(G, bisec, weight, source_modifier, sink_modifier, original_modifier, internal_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    else
        [flow, cut, reciprocalCut, matching] = Pairing(G, bisec, weight, source_modifier, sink_modifier, original_modifier); % DO FLOW, SHOULD OUTPUT SMALL SIZE OF CUT
    end
   
   % MATCHING SCALING
   matching = matching/double(match_num);
   matchrat = double(match_num)/double(match_den);
   if(~issparse(matching))
        fprintf('Matching is not sparse...\n');
   end
   fprintf(2, 'Size of cut=%d\n', length(cut));
else
   matching = sparse(1:double(n), 1:double(n), 0);
   matchrat = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT PREPARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flownumber = counter+1;

% ALL OUTPUTS ARE READY

end

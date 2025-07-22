function [A, B] = directedRoundCut(S, T, v, weights)
%DIRECTEDROUNDCUT From an initial cut, produce a nicely separate directed cut
%   

v_squared_norm = sum(v .^ 2, 2);
r_squared = median(v_squared_norm(S), Weights=double(weights(S)));
T_plus_mask = v_squared_norm(T) <= r_squared;
T_minus_mask = v_squared_norm(T) >= r_squared;
T_plus_mu = sum(weights(T(T_plus_mask)));
T_minus_mu = sum(weights(T(T_minus_mask)));
if T_plus_mu <= T_minus_mu
    A = S(v_squared_norm(S) <= r_squared);
    B = T(T_minus_mask);
else
    A = T(T_plus_mask);
    B = S(v_squared_norm(S) >= r_squared);
end


function [L, D] = laplacian(G);

D=diag(sum(G));
L=D - G;
end
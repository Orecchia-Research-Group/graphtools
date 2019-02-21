function graphstats(graphname ,t , filename)
fid = fopen(filename, 'w');
[A,n,m] = loadeg2graph(graphname);

cum = 0;
for j=1:t
  v=rand(n,1);
  [i,x] = sort(v);
  [edges, den, ex, side] = sweepcut(A, int64(x));
  ex
  cum = cum + ex;
end

D = diag(sum(A));
[V, S] = eigs(A, D, 2, 'la');
v = V(:,2);
[i,x] = sort(v);
[edges, den, ex, side] = sweepcut(A, int64(x));


fprintf(fid, 'Average Expansion\t%f\nTrials\t%d\nSpectral Gap\t%f\nSweep Cut of Eigenvector\t%f\t%d\t%d\n', cum/t, t, 1-S(2,2), ex, edges, den);

end    
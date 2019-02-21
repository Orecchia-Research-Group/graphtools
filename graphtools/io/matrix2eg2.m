function matrix2eg2(outfile, G)


fid = fopen(outfile, 'w');

if (ndims(G) ~= 2) 
	error('Not a matrix');
end

s = size(G);
n1 = s(1);
n2 = s(2);

fprintf(fid, '%d %d %d\n', n1, n2, nnz(G));
[i,j,s] = find(G);

for k=1:nnz(G)
	fprintf(fid, '%d %d %d\n', i(k)-1, j(k)-1, s(k));
end
		
fclose(fid);

end
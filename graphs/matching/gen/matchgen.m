function matching(degree, n, seed, name)

if(mod(n,2)  ~= 0)
   error('n even..\n');
end

fid = fopen( [name '.eg2'], 'w');

rand('twister', seed);

A = sparse(n,n);


for i=1:degree
   v = rand(n,1);
   [sorted , index] = sort(v);
   j = 1;
   while (j < n)
     A(index(j), index(j+1)) = 1;
     A(index(j+1), index(j))  = 1;
     j = j+2;
   end
end

k = nnz(A);
[i,j,v] = find(A);

fprintf(fid, '%d\t%d\t%d\n', n, n, k);
for t=1:k
     fprintf(fid,'%d\t%d\t1\n', i(t)-1, j(t)-1);
end

fclose(fid);

end
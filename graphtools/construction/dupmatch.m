function dupmatch(infile1, infile2,  outfile, support, num, srand)

rand('twister', srand);

[A, n1, m1] = loadeg2graph(infile1);
[B, n2, m2] = loadeg2graph(infile2);

if(support > n1 || support > n2)
  fprintf(2, 'Error: support larger than vertex set.\m');
end

G=sparse(n1 + n2, n1 + n2);
[i1, j1, s1] = find(A);
[i2, j2, s2] = find(B);
for k=1:nnz(A)
	G(i1(k), j1(k)) = s1(k);
end
for k=1:nnz(B)
	G(n1+i2(k), n1+j2(k)) = s2(k);
end

for i=1:num
   v=rand(n1,1);
   u=rand(n2,1);
   [sorted, vindex] = sort(v);
   [sorted, uindex] = sort(u);

   for j=1:support
	G(vindex(j), n1 + uindex(j)) = 1;
	G(n1 + uindex(j), vindex(j)) = 1;
   end
end

matrix2eg2(outfile, G);

end
	
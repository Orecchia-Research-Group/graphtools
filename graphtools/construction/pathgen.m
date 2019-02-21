function pathgen(outfile, n)

G=sparse(n,n);

for i=1:n-1
   G(i, i+1) = 1;
   G(i+1, i) = 1;
end

matrix2eg2(outfile, G);

end

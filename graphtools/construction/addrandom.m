function addrandom(infile, outfile, num, srand)

rand('twister', srand);

[A, n, m] = loadeg2graph(infile);

for(i=1:num)
  v = ceil(rand(1,1) * n);
  u = ceil(rand(1,1) * n);
  if(u == v)
    continue;
  end
  A(u,v) = 1;
  A(v,u) = 1;
end

matrix2eg2(outfile, A);


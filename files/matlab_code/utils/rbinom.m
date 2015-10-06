function n = rbinom(N, p)
    v = rand(N,1)<p;
    n = sum(v);
end
function idper = make_periodic_index(id, N)
    idper = mod(id-1, N)+1;
end

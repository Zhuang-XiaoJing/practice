function Ham = HNmodel_OPBC(N,beta)
t1 = 1; tm1 = 0.5; t2 = 2; tm2 = 0; % hopping amplitute
% Initialize Hamiltonian
Ham = zeros(N);
% nearest-neighbor
for ii = 1:N-1
    Ham(ii, ii+1) = t1;
    Ham(ii+1, ii) = tm1;
end
% next-nearest neighbor
for ii = 1:N-2
    Ham(ii, ii+2) = t2;
    Ham(ii+2, ii) = tm2;
end
% off-diag
Ham(1, N) = exp(-beta)*tm1;
Ham(N, 1) = exp(-beta)*t1;
Ham(1, N-1) = exp(-beta)*tm2;
Ham(N-1, 1) = exp(-beta)*t2;
Ham(2, N) = exp(-beta)*tm2;
Ham(N, 2) = exp(-beta)*t2;

end


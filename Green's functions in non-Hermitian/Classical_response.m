%%
clc;
clear;
Er = -0.96 + 1j;            % reference energy
Blist = linspace(1,90,201);  % x-axis beta
ylist = linspace(1,90,201);  % y-axis ln(|G|)
N = 100;    % Length of chain
for x = 1:length(Blist)
    
    % Solve eigen matrix and peaks data
    beta = Blist(x);
    Ham = HNmodel_OPBC(N,beta);
    evals = eig(Ham);
    
    % Green matrix
    Identity = diag(ones(1,N));
    Gmat = inv(Identity*Er - Ham);
    
    G22 = Gmat(1:2,99:100); % G Upper-Right corner blocks 2-by-2
    ylist(x) = log(abs(det(G22))); % ln |G22|
    
end
%%
ilist = zeros(1,length(Blist));
for ii = 1:length(Blist)
    beta = Blist(ii);
    elistii = eig(HNmodel_OPBC(N,beta));
    ilist(ii) = sum(abs(1./(elistii - Er)));
end
figure;
subplot(3,1,1)
plot(Blist, ilist)
subplot(3,1,2)
plot(Blist, ylist, 'o')
%%
%Derivative of ln|G1n| vs. \beta
all = [];
for ii = 2:length(Blist)-1
    dx = Blist(2) - Blist(1);
    df = ylist(ii) - ylist(ii-1);
    all = [all, df/dx];
end
%%
subplot(3,1,3)
plot(all, '-*');









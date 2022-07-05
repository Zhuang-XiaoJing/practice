%% Numerical diag
tic;
clear;
% parpool local % Starting parallel pool
gamma = sym(4/3);
t2 = 1;
T1 = linspace(-3,3,101);
N = 40;
Eall = zeros(2*N,length(T1));
% Eall = [];
parfor ii = 1:length(T1)
    t1 = T1(ii)
    cellH = [0, t1+gamma/2; t1-gamma/2, 0];
    cellH01 = [0, 0; t2, 0];
    cellH10 = [0, t2; 0, 0];
    ham = kron(diag(ones(1,N)),cellH) + kron(diag(ones(1,N-1),1),cellH01) + kron(diag(ones(1,N-1),-1),cellH10);
    e = eig(vpa(ham));
    Eall(:,ii) = e; % little bit faster
%     Eall = [Eall,e];
end
toc;
% Shut down parallel pool
% poolobj = gcp('nocreate');
% delete(poolobj);
%% Plotting
Eall = vpa(Eall);
figure;
plot(T1,abs(Eall(:,:)),'k')

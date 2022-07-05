% ------ Haldane and 3rd nearest Haldane ---------
%%
clc;clear;
Kx = linspace(0,2*pi,101);
[~,bandnum] = size(Hamiltonian(0));
band = zeros(bandnum,length(Kx));
for ii = 1:length(Kx)
    kx = Kx(ii);
    [V,D] = eig(Hamiltonian(kx));
    band(:,ii) = sort(diag(D));
end
figure;
plot(Kx,band,'k')
axis tight
%%
% ------ Functions -------
function H = Hamiltonian(kx)
Ncelly = 20;
n = 4; % 元胞中原子个数
M = 2/3; t1 = 1; t2 = 1/3; phi = pi/4;
s0 = eye(2);
sz = [1,0;0,-1];
sx = [0,1;1,0];
h00 = t1*diag(ones(1,n-1),1) + (t1*diag(ones(1,n-1),1))'+...
    t2*exp(-1j*phi)*diag(ones(1,n-2),2)+...
    (t2*exp(-1j*phi)*diag(ones(1,n-2),2))'+...
    M*kron(s0,sz);
h0y = zeros(n);
h0y(3,1) = t2*exp(1j*phi);
h0y(4,2) = t2*exp(1j*phi);
h0y(4,1) = t1;
H00 = kron(eye(Ncelly),h00) + kron(diag(ones(1,Ncelly-1),1),h0y)+...
    kron(diag(ones(1,Ncelly-1),1),h0y)'; % a slice
Hx01 = zeros(n);
Hx01(1,1) = t2*exp(1j*phi)*exp(1j*kx);
Hx01(3,3) = Hx01(1,1);
Hx01(2,2) = t2*exp(-1j*phi)*exp(1j*kx);
Hx01(4,4) = Hx01(2,2);
Hx01(2,1) = t1*exp(1j*kx);
Hx01(3,1) = t2*exp(-1j*phi)*exp(1j*kx);
Hx01(2,4) = t2*exp(1j*phi)*exp(1j*kx);
Hx01(3,4) = t1*exp(1j*kx);
Hx10 = Hx01';
Hxy01 = zeros(n);
Hxy01(3,1) = t2*exp(-1j*phi)*exp(1j*kx);
Hyx01 = zeros(n);
Hyx01(2,4) = t2*exp(1j*phi)*exp(1j*kx);
H01 = kron(eye(Ncelly),Hx01) + kron(diag(ones(1,Ncelly-1),1),Hxy01) + ...
    kron(diag(ones(1,Ncelly-1),-1),Hyx01); % hopping 0 to 1
Hxy10 = zeros(n);
Hxy10(4,2) = t2*exp(-1j*phi)*exp(-1j*kx);
Hyx10 = zeros(n);
Hyx10(1,3) = t2*exp(1j*phi)*exp(-1j*kx);
H10 = kron(eye(Ncelly),Hx10) + kron(diag(ones(1,Ncelly-1),1),Hxy10) + ...
    kron(diag(ones(1,Ncelly-1),-1),Hyx10); % hopping 1 to 0
H = H00 + H01 + H10;
end
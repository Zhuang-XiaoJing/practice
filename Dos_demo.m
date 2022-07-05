%%
% ----- Graphene ------
%%
% ----- Hamiltonian -----
clc;
clear;
t = 1;
Ncellx = 8;
Ncelly = 8;
h00 = diag(t*ones(1,3),1) + diag(t*ones(1,3),-1);
h0y = zeros(4);
h0y(4,1) = t;
H00 = kron(diag(ones(1,Ncelly)),h00) + kron(diag(ones(1,Ncelly-1),1),h0y) +...
      kron(diag(ones(1,Ncelly-1),1),h0y)';
H01 = zeros(4);
H01(2,1) = t;
H01(3,4) = t;
Hx = kron(eye(Ncelly),H01);
H = kron(eye(Ncellx),H00) + kron(diag(ones(1,Ncellx-1),1),Hx) +...
    kron(diag(ones(1,Ncellx-1),1),Hx)';
%%
[V,D] = eig(H);
%%
plot(diag(D),'o')
%%
eta = 0.05*1j;
Ef_list = linspace(-5*t,5*t,1000);
Dos_list = zeros(1,length(Ef_list));

for ii = 1:length(Ef_list)
    Ef = Ef_list(ii);
    G = inv((Ef+eta)*eye(Ncellx*Ncelly*4)-H);
    Dos = -trace(imag(G))/pi;
    Dos_list(ii) = Dos;
end
figure;
plot(Ef_list,Dos_list);



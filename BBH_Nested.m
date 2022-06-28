%  ----- BBH wilson loop and nested wilson loop -----
%%
%  -----       vx(ky)          -----
clc;
clear;
Ky = linspace(-pi,pi,20);
Kx = linspace(-pi,pi,20);
all = [];
for ii = 1:length(Ky)
    all = [all,vx(Kx,Ky(ii))];
end
figure;
plot(Ky,all,'-*')
%%
% -----     Nested py_vx       -------
clc;
clear;
Kxlist = linspace(-pi,pi,30);
Kylist = linspace(-pi,pi,30);
all = [];
for kkx = 1:length(Kxlist)
kx = Kxlist(kkx);
W = 1;
for kky = 1:length(Kylist)-1
    [V0,D0] = eig(BBH_Bulk(kx,Kylist(kky)));
    [~,ind] = sort(diag(D0));
    V0 = V0(:,ind);
    [W_occu01,W_occu02] = wannierband(Kxlist,Kylist(kky)); %对于一个ky的两条wannier band
    % use Eq.(6.5) from Benalcazar PRB, 2017
    wbasis0 = W_occu02(1)*V0(:,1) + W_occu02(2)*V0(:,2); % 这里用任一wannier band结果没有什么区别
    [V1,D1] = eig(BBH_Bulk(kx,Kylist(kky+1)));
    [~,ind] = sort(diag(D1));
    V1 = V1(:,ind);
    [W_occu11,W_occu12] = wannierband(Kxlist,Kylist(kky+1));
    % use Eq.(6.5) from Benalcazar PRB, 2017
    wbasis1 = W_occu12(1)*V1(:,1) + W_occu12(2)*V1(:,2);
    F = wbasis1'*wbasis0;
    W = F*W;
end
[Ved,Ded] = eig(BBH_Bulk(kx,Kylist(end)));
[~,ind] = sort(diag(Ded));
Ved = Ved(:,ind);
[W_occu_end1,W_occu_end2] = wannierband(Kxlist,Kylist(end));
% use Eq.(6.5) from Benalcazar PRB, 2017
wbasis_end = W_occu_end2(1)*Ved(:,1) + W_occu_end2(2)*Ved(:,2);
[Vini,Dini] = eig(BBH_Bulk(kx,Kylist(1)));
[~,ind] = sort(diag(Dini));
Vini = Vini(:,ind);
[W_occu_ini1,W_occu_ini2] = wannierband(Kxlist,Kylist(1));
% use Eq.(6.5) from Benalcazar PRB, 2017
wbasis_ini = W_occu_ini2(1)*Vini(:,1) + W_occu_ini2(2)*Vini(:,2);
Fend = wbasis_ini'*wbasis_end;
W = Fend*W;
py_vx = log(W)/(2*pi*1j);
if real(py_vx)<0
    py_vx = py_vx+1;
end
all = [all,real(py_vx)];
end
figure;
plot(Kxlist,all,'*')
ylim([0,1]);
%%
% -----------  Functions  ------------
% hamiltonian
function Ham = BBH_Bulk(kx,ky)

gamx = 0.5;
gamy = 0.5;
lamx = 1.0;
lamy = 1.0;
Ham = zeros(4);
% symmetries to quantize the quadrupole moments py_vx
x_symmetry_breaking_1 = 0.00000000;
x_symmetry_breaking_2 = 1.000000001;
y_symmetry_breaking_1 = 0.00000000;
y_symmetry_breaking_2 = 1.00000000;
Ham(1,1) = x_symmetry_breaking_1;
Ham(2,2) = y_symmetry_breaking_1;
Ham(3,3) = y_symmetry_breaking_1;
Ham(4,4) = x_symmetry_breaking_1;

Ham(1,3) = (gamx+lamx*exp(1j*kx))*y_symmetry_breaking_2;
Ham(2,4) = gamx+lamx*exp(-1j*kx);
Ham(1,4) = gamy+lamy*exp(1j*ky);
Ham(2,3) = (-gamy-lamy*exp(-1j*ky))*x_symmetry_breaking_2;
Ham(3,1) = conj(Ham(1,3));
Ham(4,2) = conj(Ham(2,4));
Ham(4,1) = conj(Ham(1,4));
Ham(3,2) = conj(Ham(2,3));
end

% 在给定ky下wannier band: vx(ky)
function ang = vx(Kxlist,ky)
Noccu = sum(eig(BBH_Bulk(0,0))<0);
W = eye(Noccu);
for jj = 1:length(Kxlist)-1
    [vec1,~] = eig(BBH_Bulk(Kxlist(jj),ky));
    [vec2,~] = eig(BBH_Bulk(Kxlist(jj+1),ky));
    F = zeros(Noccu);
    F(1,1) = vec1(:,1)'*vec2(:,1);
    F(1,2) = vec1(:,1)'*vec2(:,2);
    F(2,1) = vec1(:,2)'*vec2(:,1);
    F(2,2) = vec1(:,2)'*vec2(:,2);
    W = W*F;
end
[vecEnd,~] = eig(BBH_Bulk(Kxlist(end),ky));
[vecIni,~] = eig(BBH_Bulk(Kxlist(1),ky));
F_last = zeros(Noccu);
for i0 = 1:Noccu
    for j0 = 1:Noccu
        F_last(i0,j0) = vecEnd(:,i0)'*vecIni(:,j0);
    end
end
W = W*F_last;
ang = log(eig(W))/2/pi/1j;
for i0 = 1:2
    if real(ang(i0))<0
        ang(i0) = ang(i0)+1;
    end
end
ang = sort(real(ang));
end

% 在给定ky下vx(ky)的本征矢量wx1,wx2
function [wx1,wx2] = wannierband(Kxlist,ky)
Noccu = sum(eig(BBH_Bulk(0,0))<0);
vec1_array = [];
vec2_array = [];
for ii = 1:length(Kxlist)
    kx = Kxlist(ii);
    [eigenvector,eigenvalue] = eig(BBH_Bulk(kx,ky));
    if kx ~= pi
        vec1_array = [vec1_array,eigenvector(:,1)];
        vec2_array = [vec2_array,eigenvector(:,2)];
    else
        vec1_array = [vec1_array,vec1_array(:,1)];
        vec2_array = [vec2_array,vec2_array(:,1)];
    end
end
W_x_k = eye(Noccu);
for ii = 1:length(Kxlist)-1
    F = zeros(2);
    F(1,1) = vec1_array(:,ii+1)'*vec1_array(:,ii);
    F(2,2) = vec2_array(:,ii+1)'*vec2_array(:,ii);
    F(1,2) = vec1_array(:,ii+1)'*vec2_array(:,ii);
    F(2,1) = vec2_array(:,ii+1)'*vec1_array(:,ii);
    W_x_k = F*W_x_k;
end
[eigvec,eigval] = eig(W_x_k);
v_x = log(diag(eigval))/(2*pi*1j);
[~,ind] = sort(real(v_x));
eigvec = eigvec(:,ind);
wx1 = eigvec(:,1);
wx2 = eigvec(:,2);
end


%%
% ----- Finite system spectrum -------
clc;clear;
Ncellx = 20;
Ncelly = 20;
v = 1;
w1 = 2; w2 = 0.5;
[V,D] = eig(ham(Ncellx,Ncelly,w1,w2,v));
E_list = diag(D);
figure;
plot(E_list,'o');

%%
% ----- plot plaquette -----
% corner = abs(V(:,4*Ncellx*Ncelly/2)+V(:,4*Ncellx*Ncelly/2-1)+V(:,4*Ncellx*Ncelly/2+1)+...
%     V(:,4*Ncellx*Ncelly/2+2));
corner = 0;
for ii = 793:808
    corner = corner + abs(V(:,ii));
end

corner = reshape(corner,4,[]);
a = 10;
xlist = 0:a:(Ncellx-1)*a;
ylist = 0:-a:-(Ncelly-1)*a;
sz_center = 10;
sz = 50;
cell_position = [];
position1 = [];
sz1 = [];
position2 = [];
sz2 = [];
position3 = [];
sz3 = [];
position4 = [];
sz4 = [];
for ii = 1:length(xlist)
    x = xlist(ii);
    for jj = 1:length(ylist)
        y = ylist(jj);
        cell_position = [cell_position,[x;y]];
        %1
        position1 = [position1,[x-a/4;y+a/4]];
        sz1 = [sz1,corner(1,Ncellx*(ii-1)+jj)];
        %2
        position2 = [position2,[x+a/4;y-a/4]];
        sz2 = [sz2,corner(2,Ncellx*(ii-1)+jj)];
        %3
        position3 = [position3,[x+a/4;y+a/4]];
        sz3 = [sz3,corner(3,Ncellx*(ii-1)+jj)];
        %4
        position4 = [position4,[x-a/4;y-a/4]];
        sz4 = [sz4,corner(4,Ncellx*(ii-1)+jj)];
    end
end
sz1(sz1==0) = 1.0e-20;
sz2(find(sz2==0)) = 1.0e-20;
sz3(find(sz3==0)) = 1.0e-20;
sz4(find(sz4==0)) = 1.0e-20;
figure;
scatter(cell_position(1,:),cell_position(2,:),sz_center,'k','filled');
hold on
scatter(position1(1,:),position1(2,:),sz1*100,'b','filled');
hold on
scatter(position2(1,:),position2(2,:),sz2*100,'b','filled');
hold on
scatter(position3(1,:),position3(2,:),sz3*100,'r','filled');
hold on
scatter(position4(1,:),position4(2,:),sz4*100,'r','filled');
axis equal

%%
% ------ Density of states of finite system ---------
clc;
eta = 0.05*1j;
delta_E = 0.02;
Ef_list = -6:delta_E:6;
Dos_list = [];
parpool local
parfor ii = 1:length(Ef_list)
    ii
    Ef = Ef_list(ii);
    Gmat = inv((Ef+eta)*eye(Ncellx*Ncelly*4) - ham(Ncellx,Ncelly,w1,w2,v));
    Dos = -trace(imag(Gmat))/pi;
    Dos_list = [Dos_list,Dos];
end
% Shut down parallel pool
poolobj = gcp('nocreate');
delete(poolobj);

figure;
plot(Ef_list,Dos_list);

%%
% ---- 3D Bulk Band ---------
clc;clear;
v = 1;
w1 = 2; w2 = 0;
num = 50;
band_num = size(hamBulk(0,0,w1,w2,v),1);
Kx = linspace(-pi,pi,50);
Ky = linspace(-pi,pi,50);
band = zeros(num,num,band_num);
for ii = 1:length(Kx)
    kx = Kx(ii);
    for jj = 1:length(Ky)
        ky = Ky(jj);
        [V,D] = eig(hamBulk(kx,ky,w1,w2,v));
        band(ii,jj,:) = sort(diag(D));
    end
end
figure;
for index = 1:band_num
    mesh(band(:,:,index));
    hold on
end
%%
% ------ High symmetry line G-Corner-
clc;clear;
v=1; w1 = 2; w2 = 0;
step = 0.001;
b = 2*pi;
path1 = 0:step:sqrt(2);
path2 = (sqrt(2)+step):step:(1+sqrt(2));
path3 = (1+sqrt(2)+step):step:(2+sqrt(2));
E_array = [];
for ii = 1:length(path1)
    kx = path1(ii)*(b/2)/sqrt(2);
    ky = path1(ii)*(b/2)/sqrt(2);
    [V,D] = eig(hamBulk(kx,ky,w1,w2,v));
    E_array = [E_array,sort(diag(D))];
end
for ii = 1:length(path2)
    kx = (1+sqrt(2)-path2(ii))*(b/2);
    ky = b/2;
    [V,D] = eig(hamBulk(kx,ky,w1,w2,v));
    E_array = [E_array,sort(diag(D))];
end
for ii = 1:length(path3)
    kx = 0;
    ky = (b/2)*(2+sqrt(2)-path3(ii));
    [V,D] = eig(hamBulk(kx,ky,w1,w2,v));
    E_array = [E_array,sort(diag(D))];
end
figure;
plot(cat(2,path1,path2,path3),E_array,'k')
%%
% ------ Functions ---------
% ------ Finite Hamiltonian -----
function H = ham(Ncellx,Ncelly,w1,w2,v)
vx = v; vy = v; w1x = w1; w1y = w1; w2x = w2; w2y = w2;
h00 = zeros(4);
h00(1:2,3:4) = [-vx,vy;vy,vx];
h00(3:4,1:2) = h00(1:2,3:4)';
h01 = zeros(4);
h01(2,4) = w1x;
h01(3,1) = -w1x;
h02 = zeros(4);
h02(2,4) = w2x;
h02(3,1) = -w2x;
M00 = kron(eye(Ncellx),h00) + kron(diag(ones(1,Ncellx-1),1),h01)' +...
    kron(diag(ones(1,Ncellx-2),2),h02)' + kron(diag(ones(1,Ncellx-1),1),h01) +...
    kron(diag(ones(1,Ncellx-2),2),h02);
hy01 = zeros(4);
hy01(2,3) = w1y;
hy01(4,1) = w1y;
M01 = kron(eye(Ncellx),hy01);
hy02 = zeros(4);
hy02(2,3) = w2y;
hy02(4,1) = w2y;
M02 = kron(eye(Ncellx),hy02);
H = kron(eye(Ncelly),M00) + kron(diag(ones(1,Ncelly-1),1),M01) +...
    kron(diag(ones(1,Ncelly-1),1),M01)' + kron(diag(ones(1,Ncelly-2),2),M02) + ...
    kron(diag(ones(1,Ncelly-2),2),M02)';
end

% ------ Bulk Hamiltonian -----

function H = hamBulk(kx,ky,w1,w2,v)
vx = v; vy = v; w1x = w1; w1y = w1; w2x = w2; w2y = w2;
hQTI = zeros(2);
hQTI(1,1) = -vx - w1x*exp(-1j*kx);
hQTI(1,2) = vy + w1y*exp(1j*ky);
hQTI(2,1) = vy + w1y*exp(-1j*ky);
hQTI(2,2) = vx + w1x*exp(1j*kx);
hSLR = zeros(2);
hSLR(1,1) = -w2x*exp(-2*1j*kx);
hSLR(1,2) = w2y*exp(2*1j*ky);
hSLR(2,1) = w2y*exp(-2*1j*ky);
hSLR(2,2) = w2x*exp(2*1j*kx);
H = kron([0,1;0,0],hQTI+hSLR) + kron([0,1;0,0],hQTI+hSLR)';
end

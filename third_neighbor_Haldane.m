% ------ Haldane and 3rd nearest Haldane ---------
%%
% ------- strip band ----------
clc;clear;
Ncelly = 40;
Kx = linspace(0,2*pi,101);
[~,bandnum] = size(Hamiltonian(0,Ncelly));
band = zeros(bandnum,length(Kx));
for ii = 1:length(Kx)
    kx = Kx(ii);
    [V,D] = eig(Hamiltonian(kx,Ncelly));
    band(:,ii) = sort(diag(D));
end
figure;
plot(Kx,band,'k',Kx,band(Ncelly*2,:),'r',Kx,band(Ncelly*2+1,:),'r')
axis tight

% figure;
% kx = pi+0.1;
% [V,D] = eig(Hamiltonian(kx,Ncelly));
% plot(abs(V(:,Ncelly*2+1)))
%%
% ---- bulk band -------
clc;clear;
M = 0;
phi = pi/2;
Kx = linspace(-pi,pi,101);
Ky = linspace(-pi,pi,101);
[~,band_num] = size(haldaneBulk(0,0,M,phi));
band = zeros(length(Kx),length(Ky),band_num);
for ii = 1:length(Kx)
    kx = Kx(ii);
    for jj = 1:length(Ky)
        ky = Ky(jj);
        [V,D] = eig(haldaneBulk(kx,ky,M,phi));
        band(ii,jj,:) = sort(real(diag(D)));
    end
end
figure;
[X,Y] = meshgrid(Kx,Ky);
for ind = 1:band_num
    mesh(X,Y,band(:,:,ind));
    hold on
end
%%
% ---- phase diagram ----
clc;clear;
num = 50;
a = sqrt(3);
G = 4*pi/sqrt(3)/a;
M_list = linspace(-6,6,51);
phi_list = linspace(-pi,pi,51);
K1 = linspace(0,1,num); % steps in 1st Brillouin zone
K2 = linspace(0,1,num);
Ch_all = zeros(length(M_list),length(phi_list)); % initialize Chern # matrix
% Create waitbar
h = waitbar(0,'Please wait ... ');
steps1 = length(M_list);
steps2 = length(phi_list);
steps = steps1*steps2;
% calculation begins:
for m = 1:steps1
    for p = 1:steps2
        M = M_list(m);
        phi = phi_list(p);
        F_all = zeros(length(K1),length(K2));
        for ii = 1:length(K1)-1
            for jj = 1:length(K2)-1
                % transfer to diamond 1st BZ
                kx = sqrt(3)/2*G*K1(ii)+sqrt(3)/2*G*K2(jj);   
                ky = -G/2*K1(ii) + G/2*K2(jj);
                kxx = sqrt(3)/2*G*K1(ii+1)+sqrt(3)/2*G*K2(jj);
                kyy = -G/2*K1(ii) + G/2*K2(jj+1);
                % Calculate connections
                V = getvec(haldaneBulk(kx,ky,M,phi));
                Vkx = getvec(haldaneBulk(kxx, ky,M,phi));
                Vky = getvec(haldaneBulk(kx,kyy,M,phi));
                Vkxky  = getvec(haldaneBulk(kxx,kyy,M,phi));
                U1 = V'*Vkx/abs(V'*Vkx);
                U12 = Vkx'*Vkxky/abs(Vkx'*Vkxky);
                U21 = Vky'*Vkxky/abs(Vky'*Vkxky);
                U2 = V'*Vky/abs(V'*Vky);
                % Japanese way
                F = log(U1*U12/U21/U2);
                F_all(ii,jj) = F;
            end
        end
        % summerize Berry curvature F at one point on the phase diagram to
        % obtain Chern #
        Ch_all(m,p) = sum(sum(F_all)); 
    end
    % current process
    per  = m*p / steps;
    waitbar(per,h,sprintf('%2.0f%%',per*100));
end
close(h)
% plot
figure;
imagesc(imag(Ch_all)/3.14);
ax = gca;
ax.YDir = 'normal';
colorbar;
%%
% ----- Bulk Hamiltonian -------
function h = haldaneBulk(kx,ky,M,phi)

b1 = [-sqrt(3)/2,3/2];
b2 = [-sqrt(3)/2,-3/2];
b3 = [sqrt(3),0];
s0 = eye(2); sx = [0,1;1,0]; sy = [0,-1j;1j,0]; sz = [1,0;0,-1];
k = [kx,ky];
t1 = 1;
t2 = 1;
t3 = 5;
h0 = 2*t2*cos(phi)*(cos(dot(k,b1))+cos(dot(k,b2))+cos(dot(k,b3)));
hx = t1*(1+cos(dot(k,b1))+cos(dot(k,b2))) + t3*(2*cos(dot(k,b1+b2))+cos(dot(k,b1-b2)));
hy = t1*(sin(dot(k,b1))-sin(dot(k,b2))) + t3*sin(dot(k,b1-b2));
hz = M - 2*t2*sin(phi)*(sin(dot(k,b1))+sin(dot(k,b2))+sin(dot(k,b3)));
h = h0*s0 + hx*sx + hy*sy + hz*sz;

end

% ----- Get occupied band ----
function occu = getvec(H)
[V,D] = eig(H);
[~,ind] = sort(diag(D),'ascend');
occu = V(:,ind(1)); % only one occupied band
end
%%
% ------ strip hamiltonian -------
function H = Hamiltonian(kx,Ncelly)

n = 4; % 元胞中原子个数
M = 0; t1 = 1; t2 = 1/3; phi = pi/2; t3 = 0;
s0 = eye(2);
sz = [1,0;0,-1];
sx = [0,1;1,0];
h00 = t1*diag(ones(1,n-1),1) + (t1*diag(ones(1,n-1),1))'+...
    t2*exp(-1j*phi)*diag(ones(1,n-2),2)+...
    (t2*exp(-1j*phi)*diag(ones(1,n-2),2))'+...
    M*kron(s0,sz)+...
    [0,0,0,t3;0,0,0,0;0,0,0,0;t3',0,0,0];
h0y = zeros(n);
h0y(3,1) = t2*exp(1j*phi);
h0y(3,2) = t3;
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
Hx01(2,3) = t3*exp(1j*kx);
Hx01(3,1) = t2*exp(-1j*phi)*exp(1j*kx);
Hx01(3,2) = t3*exp(1j*kx);
Hx01(2,4) = t2*exp(1j*phi)*exp(1j*kx);
Hx01(3,4) = t1*exp(1j*kx);
Hx10 = Hx01';
Hxy01 = zeros(n);
Hxy01(3,1) = t2*exp(-1j*phi)*exp(1j*kx);
Hxy01(4,1) = t3*exp(1j*kx);
Hyx01 = zeros(n);
Hyx01(2,4) = t2*exp(1j*phi)*exp(1j*kx);
Hyx01(1,4) = t3*exp(1j*kx);
H01 = kron(eye(Ncelly),Hx01) + kron(diag(ones(1,Ncelly-1),1),Hxy01) + ...
    kron(diag(ones(1,Ncelly-1),-1),Hyx01); % hopping 0 to 1
Hxy10 = zeros(n);
Hxy10(4,2) = t2*exp(-1j*phi)*exp(-1j*kx);
Hxy10(4,1) = t3*exp(-1j*kx);
Hyx10 = zeros(n);
Hyx10(1,3) = t2*exp(1j*phi)*exp(-1j*kx);
Hyx10(1,4) = t3*exp(-1j*kx);
H10 = kron(eye(Ncelly),Hx10) + kron(diag(ones(1,Ncelly-1),1),Hxy10) + ...
    kron(diag(ones(1,Ncelly-1),-1),Hyx10); % hopping 1 to 0
H = H00 + H01 + H10;
end
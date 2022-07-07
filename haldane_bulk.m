% --- Haldane bulk ---
clc;clear;
Kx = linspace(-pi,pi,101);
Ky = linspace(-pi,pi,101);
[~,band_num] = size(haldaneBulk(0,0));
band = zeros(length(Kx),length(Ky),band_num);
for ii = 1:length(Kx)
    kx = Kx(ii);
    for jj = 1:length(Ky)
        ky = Ky(jj);
        [V,D] = eig(haldaneBulk(kx,ky));
        band(ii,jj,:) = sort(diag(D));
    end
end
figure;
[X,Y] = meshgrid(Kx,Ky);
for ind = 1:band_num
    mesh(X,Y,band(:,:,ind));
    hold on
end
%%
% --- Berry curvature ---
clc;clear;
num = 50;
a = sqrt(3);
G = 4*pi/sqrt(3)/a;
M = 0.5;
phi = -pi/2;
K1 = linspace(0,1,num);
K2 = linspace(0,1,num);

F_all = zeros(length(K1),length(K2));
for ii = 1:length(K1)-1
    for jj = 1:length(K2)-1
        kx = sqrt(3)/2*G*K1(ii)+sqrt(3)/2*G*K2(jj);
        ky = -G/2*K1(ii) + G/2*K2(jj);
        kxx = sqrt(3)/2*G*K1(ii+1)+sqrt(3)/2*G*K2(jj);
        kyy = -G/2*K1(ii) + G/2*K2(jj+1);
        V = getvec(haldaneBulk(kx,ky,M,phi));
        Vkx = getvec(haldaneBulk(kxx, ky,M,phi));
        Vky = getvec(haldaneBulk(kx,kyy,M,phi));
        Vkxky  = getvec(haldaneBulk(kxx,kyy,M,phi));
        U1 = V'*Vkx/abs(V'*Vkx);
        U12 = Vkx'*Vkxky/abs(Vkx'*Vkxky);
        U21 = Vky'*Vkxky/abs(Vky'*Vkxky);
        U2 = V'*Vky/abs(V'*Vky);
        F = log(U1*U12/U21/U2);
        F_all(ii,jj) = F;
    end
end

figure;
imagesc(imag(F_all)+real(F_all));
ax = gca;
ax.YDir = 'normal';
shading interp
colorbar;
Ch = (sum(sum(F_all)))
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
% ----- Hamiltonian -------
function h = haldaneBulk(kx,ky,M,phi)

b1 = [-sqrt(3)/2,3/2];
b2 = [-sqrt(3)/2,-3/2];
b3 = [sqrt(3),0];
k = [kx,ky];
s0 = eye(2); sx = [0,1;1,0]; sy = [0,-1j;1j,0]; sz = [1,0;0,-1];
% M = 0.5;
t1 = -1;
t2 = -1;
% phi = -pi/2;
h0 = 2*t2*cos(phi)*(cos(dot(k,b1))+cos(dot(k,b2))+cos(dot(k,b3)));
hx = t1*(1 + cos(dot(k,b1)) + cos(dot(k,b2)));
hy = t1*(sin(dot(k,b1)) - sin(dot(k,b2)));
hz = M - 2*t2*sin(phi)*(sin(dot(k,b1)) + sin(dot(k,b2)) + sin(dot(k,b3)));
h = h0*s0 + hx*sx + hy*sy + hz*sz;

end
% ----- Get occupied band ----
function occu = getvec(H)
[V,D] = eig(H);
[~,ind] = sort(diag(D),'ascend');
occu = V(:,ind(1)); % only one occupied band
end
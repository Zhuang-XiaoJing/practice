%% bulk
clear;
a = 1;  % lattice constant
b = 2*pi/a; % Bloch vector G
% define the BZ path
path1 = 0:0.01:1;
path2 = 1:0.01:2;
path3 = 2:0.01:2+sqrt(2);
eigvals = [];   % initialize the band matrix
% Gamma---X
for ii=1:length(path1)
    kx = path1(ii)*b/2;
    ky = 0;
    H = hamD4(kx,ky);
    eigvals = [eigvals,sort(eig(H))];
end
% X---M
for ii = 1:length(path2)
    kx = b/2;
    ky = (path2(ii)-1)*b/2;
    H = hamD4(kx,ky);
    eigvals = [eigvals, sort(eig(H))];
end
% M---Gamma
for ii = 1:length(path3)
    kx = b/2*(2+sqrt(2)-path3(ii))/sqrt(2);
    ky = kx;
    H = hamD4(kx,ky);
    eigvals = [eigvals, sort(eig(H))];
end
% Join all paths and make a plot
path = cat(2,path1,path2,path3);
figure
plot(path,eigvals)

%% nano square disk with size 5*5 sites
clear;clc;
% define the side length of disk
N = 9;  % odd number
% define the hopping
J = -1; K = -4; M = -4;

% in-cell matrix
off_diag00 = kron(ones(1,(N-1)/2),[K,J]);
h00 = diag(off_diag00,1) + diag(off_diag00,-1);

% KM hopping matrix
off_diagKM = kron(ones(1,(N-1)/2),[M,0]);
hKM = diag(off_diagKM,1) + diag(off_diagKM,-1) + diag(K*ones(1,N));

% J hopping matrix
hJ = diag(J*ones(1,N));

% Hamiltonian all
pos_hKM = kron(ones(1,(N-1)/2),[1,0]);
pos_hJ = kron(ones(1,(N-1)/2),[0,1]);
pos_h00 = ones(1,N);
Hdisk_off = kron(diag(pos_hKM,1),hKM) + kron(diag(pos_hJ,1),hJ);
Hdisk_diag = kron(diag(pos_h00),h00);
Hdisk = Hdisk_off + Hdisk_off' + Hdisk_diag;

% visualize matrix
figure;
imagesc(Hdisk);
colorbar
axis equal
axis tight
% set(gca,'ydir','normal');

[eigenvec,eigenval] = eig(Hdisk);

% eigen scatter plot
figure;
index = 1:N*N;
elist = sort(diag(eigenval));
scatter(index,elist,'filled')

% visualize corner state
CornerState = eigenvec(:,25).^2;
corner = reshape(CornerState,[N,N]);
figure;
imagesc(corner);
colorbar;
set(gca,'ydir','normal');
axis equal
axis tight
% corner = zeros(N,N);
% for y=N:-1:1
%     flag = 1;
%     for x=1:N
%         corner(x,y) = CornerState(flag);
%     end
%     flag = flag + 1;
% end

% visualize edge states
EdgeState = eigenvec(:,23).^2;
edge = reshape(EdgeState,[N,N]);
figure;
imagesc(edge);
colorbar;
set(gca,'ydir','normal');
axis equal
axis tight
%%
% Spectrum of nano disk of 9*9 structure with J=1 and K=4 changing M
clear;clc;close;
% 固定参数
J = -1; K = -4; N = 9;
% 变化参数
MJlist = linspace(0,8,101);
% 初始化矩阵
emat = [];
cmat = zeros(length(MJlist), N^2);
tic;
% 储存能带
for ii = 1: length(MJlist)
    M = MJlist(ii)*J;
    [eigenvec,eigenval] = eig(hamD4disk(M));
    emat = [emat,real(diag(eigenval))];
    for jj = 1: length(diag(eigenval))
        cmat(ii,jj) = sum(eigenvec(:,jj).^4);
    end
end

% 画图
figure;
for ii = 1:N^2
    Ei = emat(ii,:); % y = sin(x)
    Ei(end) = NaN;   % y(end) = NaN;
    ci = cmat(:,ii);
    ci(end) = NaN;
    c = ci;
    % patch(x, y, c, 'EdgeColor'...)
    patch(MJlist,Ei,c, 'EdgeColor', 'interp', 'LineWidth', 3);
%     scatter(MJlist,Ei,'o');
    hold on
end
colorbar
toc;
% Visualize edge states
%%
figure;
M = -5.085;
[eigenvec,eigenval] = eig(hamD4disk(M));
imagesc(reshape(eigenvec(:,68).^2,[N,N]));
set(gca,'ydir','normal');
axis equal
axis tight
axis square
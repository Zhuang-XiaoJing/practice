%%
%-------------------Demo of sorting eigsystems----------------------------
clear;
A = magic(5);
[V,D] = eig(A);
[d, ind] = sort(diag(D))
Ds = D(ind,ind)
Vs = V(:,ind)
%%
%----------------------General way by YuRui--------------------------------
clc;clear;
Ky = linspace(0, 2*pi, 500);
Kx = linspace(-pi, pi, 50); %每条wilson loop上的点数
u = 1; g = 0.1; C = 0.1*[0,-1j;1j,0];
band_num = sum(eig(BHZ2(0,0))<0); %计算占据态的个数
all = [];
for y0 = 1:length(Ky)
    wmat = eye(band_num);
    for x0 = 1:length(Kx)-1
        F = zeros(band_num);
%         [V1,D1] = eig(BHZ2(Kx(x0),Ky(y0)));   % kx处的本征
        [V1,D1] = eig(BHZ(Kx(x0),Ky(y0),u,C,g));   % kx处的本征
        [d1,ind1] = sort(diag(D1));
        Vs1 = V1(:,ind1);
%         [V2,D2] = eig(BHZ2(Kx(x0+1),Ky(y0))); % kx+dx处的本征
        [V2,D2] = eig(BHZ(Kx(x0+1),Ky(y0),u,C,g));   % kx+dx处的本征
        [d2,ind2] = sort(diag(D2));
        Vs2 = V2(:,ind2);
        for ii = 1:band_num
            for jj=1:band_num
                F(ii,jj) = Vs1(:,ii)'*Vs2(:,jj);
            end
        end
        wmat = wmat*F;
    end
    Fend1 = zeros(band_num);
%     [V,D] = eig(BHZ2(Kx(end),Ky(y0)));  %尾部本征
    [V,D] = eig(BHZ(Kx(end),Ky(y0),u,C,g));  %尾部本征
    [d,ind] = sort(diag(D));
    Vs = V(:,ind);
%     [V0,D0] = eig(BHZ2(Kx(1),Ky(y0)));  %头部本征
    [V0,D0] = eig(BHZ(Kx(1),Ky(y0),u,C,g));  %头部本征
    [d0,ind0] = sort(diag(D0));
    Vs0 = V0(:,ind);
    for ii = 1:band_num
        for jj = 1:band_num
            Fend1(ii,jj) = Vs(:,ii)'*Vs0(:,jj);
        end
    end
    wmat = wmat*Fend1;
    all = [all,sort(angle(eig(wmat)))]; % 这里注意是angle不是imag！！
end

figure;
plot(Ky, all, 'k*')
xlim([Ky(1),Ky(end)])
ylim([-pi,pi])
xlabel('$k_y$','interpreter','latex','fontsize',12)
ylabel('Wannier Center $\theta$','Interpreter','latex','fontsize',12)
grid on

%%
%-----------------General way by YuRui for Graphene------------------------
clc;clear;
Ky = linspace(0, 2*pi, 500);
Kx = linspace(-pi, pi, 50); %每条wilson loop上的点数
u = 1; g = 0.1; C = 0.1*[0,-1j;1j,0];
band_num = sum(eig(KMBulk(0,0))<0); %计算占据态的个数
all = [];
for y0 = 1:length(Ky)
    wmat = eye(band_num);
    for x0 = 1:length(Kx)-1
        F = zeros(band_num);
        [V1,D1] = eig(KMBulk(Kx(x0),Ky(y0)));   % kx处的本征
        [d1,ind1] = sort(diag(D1));
        Vs1 = V1(:,ind1);
        [V2,D2] = eig(KMBulk(Kx(x0+1),Ky(y0))); % kx+dx处的本征
        [d2,ind2] = sort(diag(D2));
        Vs2 = V2(:,ind2);
        for ii = 1:band_num
            for jj=1:band_num
                F(ii,jj) = Vs1(:,ii)'*Vs2(:,jj);
            end
        end
        wmat = wmat*F;
    end
    Fend1 = zeros(band_num);
    [V,D] = eig(KMBulk(Kx(end),Ky(y0)));  %尾部本征
    [d,ind] = sort(diag(D));
    Vs = V(:,ind);
    [V0,D0] = eig(KMBulk(Kx(1),Ky(y0)));  %头部本征
    [d0,ind0] = sort(diag(D0));
    Vs0 = V0(:,ind);
    for ii = 1:band_num
        for jj = 1:band_num
            Fend1(ii,jj) = Vs(:,ii)'*Vs0(:,jj);
        end
    end
    wmat = wmat*Fend1;
    all = [all,sort(angle(eig(wmat)))]; % 这里注意是angle不是imag！！
end
figure;
plot(Ky, all, 'k*')
xlim([Ky(1),Ky(end)])
ylim([-pi,pi])
xlabel('$k_y$','interpreter','latex','fontsize',12)
ylabel('Wannier Center $\theta$','Interpreter','latex','fontsize',12)
grid on

%%
%---------------------The way used by László Oroszlány---------------------
clc;clear;
C = 0.1*[0,-1j;1j,0];
Ky = linspace(0, pi, 500);
Kx = linspace(-pi, pi, 10); %每条wilson loop上的点数
for y0 = 1:length(Ky)
    Wmat = eye(4); %初始化
    % 开始kky处的一条Loop
    for x0 = 1:length(Kx)
        H = BHZ(Kx(x0),Ky(y0),-2.1,C,0);
        [V,D] = eig(H);
        P = V(:,diag(D)<0)*V(:,diag(D)<0)';
        Wmat = Wmat*P; % 累乘
    end
    % 结束一条Loop
    Weig = eig(Wmat); 
    index = find(Weig>1e-10);
    plot(Ky(y0)*ones(1,length(index)),angle(Weig(index)),'b*');
    hold on
end
xlim([Ky(1),Ky(end)])
ylim([-pi,pi])
xlabel('ky')
ylabel('Wannier Center')









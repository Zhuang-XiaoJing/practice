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
Kx = linspace(-pi, pi, 50); %ÿ��wilson loop�ϵĵ���
u = 1; g = 0.1; C = 0.1*[0,-1j;1j,0];
band_num = sum(eig(BHZ2(0,0))<0); %����ռ��̬�ĸ���
all = [];
for y0 = 1:length(Ky)
    wmat = eye(band_num);
    for x0 = 1:length(Kx)-1
        F = zeros(band_num);
%         [V1,D1] = eig(BHZ2(Kx(x0),Ky(y0)));   % kx���ı���
        [V1,D1] = eig(BHZ(Kx(x0),Ky(y0),u,C,g));   % kx���ı���
        [d1,ind1] = sort(diag(D1));
        Vs1 = V1(:,ind1);
%         [V2,D2] = eig(BHZ2(Kx(x0+1),Ky(y0))); % kx+dx���ı���
        [V2,D2] = eig(BHZ(Kx(x0+1),Ky(y0),u,C,g));   % kx+dx���ı���
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
%     [V,D] = eig(BHZ2(Kx(end),Ky(y0)));  %β������
    [V,D] = eig(BHZ(Kx(end),Ky(y0),u,C,g));  %β������
    [d,ind] = sort(diag(D));
    Vs = V(:,ind);
%     [V0,D0] = eig(BHZ2(Kx(1),Ky(y0)));  %ͷ������
    [V0,D0] = eig(BHZ(Kx(1),Ky(y0),u,C,g));  %ͷ������
    [d0,ind0] = sort(diag(D0));
    Vs0 = V0(:,ind);
    for ii = 1:band_num
        for jj = 1:band_num
            Fend1(ii,jj) = Vs(:,ii)'*Vs0(:,jj);
        end
    end
    wmat = wmat*Fend1;
    all = [all,sort(angle(eig(wmat)))]; % ����ע����angle����imag����
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
Kx = linspace(-pi, pi, 50); %ÿ��wilson loop�ϵĵ���
u = 1; g = 0.1; C = 0.1*[0,-1j;1j,0];
band_num = sum(eig(KMBulk(0,0))<0); %����ռ��̬�ĸ���
all = [];
for y0 = 1:length(Ky)
    wmat = eye(band_num);
    for x0 = 1:length(Kx)-1
        F = zeros(band_num);
        [V1,D1] = eig(KMBulk(Kx(x0),Ky(y0)));   % kx���ı���
        [d1,ind1] = sort(diag(D1));
        Vs1 = V1(:,ind1);
        [V2,D2] = eig(KMBulk(Kx(x0+1),Ky(y0))); % kx+dx���ı���
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
    [V,D] = eig(KMBulk(Kx(end),Ky(y0)));  %β������
    [d,ind] = sort(diag(D));
    Vs = V(:,ind);
    [V0,D0] = eig(KMBulk(Kx(1),Ky(y0)));  %ͷ������
    [d0,ind0] = sort(diag(D0));
    Vs0 = V0(:,ind);
    for ii = 1:band_num
        for jj = 1:band_num
            Fend1(ii,jj) = Vs(:,ii)'*Vs0(:,jj);
        end
    end
    wmat = wmat*Fend1;
    all = [all,sort(angle(eig(wmat)))]; % ����ע����angle����imag����
end
figure;
plot(Ky, all, 'k*')
xlim([Ky(1),Ky(end)])
ylim([-pi,pi])
xlabel('$k_y$','interpreter','latex','fontsize',12)
ylabel('Wannier Center $\theta$','Interpreter','latex','fontsize',12)
grid on

%%
%---------------------The way used by L��szl�� Oroszl��ny---------------------
clc;clear;
C = 0.1*[0,-1j;1j,0];
Ky = linspace(0, pi, 500);
Kx = linspace(-pi, pi, 10); %ÿ��wilson loop�ϵĵ���
for y0 = 1:length(Ky)
    Wmat = eye(4); %��ʼ��
    % ��ʼkky����һ��Loop
    for x0 = 1:length(Kx)
        H = BHZ(Kx(x0),Ky(y0),-2.1,C,0);
        [V,D] = eig(H);
        P = V(:,diag(D)<0)*V(:,diag(D)<0)';
        Wmat = Wmat*P; % �۳�
    end
    % ����һ��Loop
    Weig = eig(Wmat); 
    index = find(Weig>1e-10);
    plot(Ky(y0)*ones(1,length(index)),angle(Weig(index)),'b*');
    hold on
end
xlim([Ky(1),Ky(end)])
ylim([-pi,pi])
xlabel('ky')
ylabel('Wannier Center')









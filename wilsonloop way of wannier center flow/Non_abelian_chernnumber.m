% Wilson loop相交就不可辨识，有gap也不可以辨识，需要继续处理
% 注意要单独算能带，否则在简并点算两次结果不同
% 这里算的是Wannier Centers,也称为Wannier bands
clear; clc; tic; 

e_all_real=load('Non_abelian_wannier_square_first_real.txt');
e_all_imag=load('Non_abelian_wannier_square_first_imag.txt');
e_all_real1=e_all_real(:,4:885);
e_all_imag1=e_all_imag(:,4:885);
e_all=e_all_real1+i.*e_all_imag1;

N_mash=13467;
Ny=20; dky=4/sqrt(3)/Ny; kys=-2/sqrt(3):dky:2/sqrt(3); 
Nx=Ny; dkx=2/Nx; kxs=-1:dkx:1; %四方

Estates=zeros(2*N_mash,1,length(kxs),length(kys));

for ikx=1:length(kxs)
    ikx
    for iky=1:length(kys)
        
    e1=e_all(:,42*(ikx-1)+2*iky-1);
    e2=e_all(:,42*(ikx-1)+2*iky);
    e_12=[e1;e2];
    
    Estates(:,:,ikx,iky)=e_12;
    end
end


all=[];
C_number=0;
for ikx=1:length(kxs)-1
    ikx 
    for iky=1:length(kys)-1
    A12=Estates(:,1,ikx,iky)'*Estates(:,1,ikx+1,iky); 
    A23=Estates(:,1,ikx+1,iky)'*Estates(:,1,ikx+1,iky+1);
    A34=Estates(:,1,ikx+1,iky+1)'*Estates(:,1,ikx,iky+1);
    A41=Estates(:,1,ikx,iky+1)'*Estates(:,1,ikx,iky);
    C_number=C_number+imag(log(A12*A23*A34*A41)); 
    end
    
end
C_number=-C_number/(2*pi)

time=toc



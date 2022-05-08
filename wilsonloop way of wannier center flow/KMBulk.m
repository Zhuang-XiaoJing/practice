function H = KMBulk(kx,ky)
t=-1;
gamma_so=0.06*t;
gamma_v=0.01*t;
gamma_r=0.05*t;

s0=eye(2);
sx=[0,1;1,0];
sy=[0,-1j;1j,0];
sz=[1,0;0,-1];
s1=[0,0;1,0];
s2=[0,1;0,0];

k=[kx,ky];
v1=[sqrt(3),0]; v3=[sqrt(3)/2,3/2]; v2=[sqrt(3)/2,-3/2];

H00=kron(sz,gamma_v*s0)+kron(sx,t*s0)+kron(sy,gamma_r*sx);      
H01=kron(-sz,1j*gamma_so*sz);
H1=H01*exp(1j*dot(k,v1));
H02=kron(-sz,1j*gamma_so*sz)+kron(s1,t*s0+1j*gamma_r*(-1/2*sx-sqrt(3)/2*sy));
H2=H02*exp(1j*dot(k,v2));
H03=kron(-sz,1j*gamma_so*sz)+kron(s2,t*s0+1j*gamma_r*(1/2*sx-sqrt(3)/2*sy));
H3=H03*exp(1j*dot(k,v3));
        
H=H00+H1+H1'+H2+H2'+H3+H3';
end


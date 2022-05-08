function H = BHZ(kx,ky,u,C,g)
s0 = [1,0;0,1];
s1 = [0,1;1,0];
s2 = [0,-1j;1j,0];
s3 = [1,0;0,-1];
H = kron(s0, (u+cos(kx)+cos(ky))*s3 + sin(ky)*s2)...
    + kron(s3, sin(kx)*s1)...
    + kron(s1, C)...
    + kron(g*s3, s2*(cos(kx)+cos(7*ky)-2));
end


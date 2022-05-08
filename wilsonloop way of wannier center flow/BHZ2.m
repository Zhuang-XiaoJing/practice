function H = BHZ2(kx,ky)
A = -13.68; B =-16.9; C = -0.0263; D = -0.514; M = 0.1; delta= 1.20;
sx=[0, 1; 1, 0];sy=1i*[0, -1;1, 0];sz=[1, 0;0, -1];
ek=C-2*D*(2-cos(kx)-cos(ky)); 
d1=A*sin(kx); 
d2=A*sin(ky); 
d3=-2*B*(2-M/(2*B)-cos(kx)-cos(ky));
H1=ek+d1*sx+d2*sy+d3*sz; H2=conj(ek-d1*sx-d2*sy+d3*sz);
H=[H1,zeros(2);zeros(2),H2];
end


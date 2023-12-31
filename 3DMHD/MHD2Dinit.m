% Initialize the state and run the 2D MHD solver

B0 = 20;
mi = 1.67e-24;
T0 = 4.8e-12;
cs = sqrt(2*T0/mi);
n0 = 4e8;
Va = B0/sqrt(4*pi*n0*mi);
beta = cs^2/Va^2;
j0 = 75;i0 = 5;D = 14; clear u0; clear u;
for i = 1:100
for j = 1:200
u0(i,j,1) = 2e8*n0^-1;
u0(i,j,2) = 0;
u0(i,j,3) = 0;
u0(i,j,4) = u0(i,j,1)*4.8e-12*T0^-1*beta;
u0(i,j,5) = 20.0*B0^-1*(1+i0)*(i+i0)./((i+i0).^2+(j-100).^2);
u0(i,j,6) = 20.0*B0^-1*(1+i0)*(-(j-100))./((i+i0).^2+(j-100).^2);
end
end
for j = 1:100
u0(1,j,3) = 3e7*Va^-1*(exp(-(j-round(200/2*(1/4.5+1/2.5))).^2/D^2));
u0(1,j,1) = 2e8*n0^-1*(1+(exp(-(j-round(200/2*(1/4.5+1/2.5))).^2/D^2)));
u0(1,j,4) = u0(1,j,1)*4.8e-12*T0^-1*beta;
end
tic, [u t] = MHD_MacCormack8(7000,10,1/5,5/3,u0); toc

%% initialize parameters
clear; clc; close all;
f0 = 10e9;
lambda0 = 3e8/f0;
f = linspace(0.2*f0,f0,50);
dx = lambda0/2;
dy = dx;
wd = lambda0/5;
deltad = wd;
hd = lambda0/4;
Nx = 7;
Ndum = 2;
lambda = 3e8./f;
k0 = 2*pi./lambda;
% th1 = [eps pi/4];
% ph1 = [eps pi/2-eps];
% [th, ph] = meshgrid(th1, ph1);
th=eps;ph = eps;
my_lim = 20;
v0=1;
Zl = 400;
%% find Y matrix
Y1 = zeros(Nx+2*Ndum,1);
v = zeros(Nx+2*Ndum,1);
Znx = zeros(length(f), Nx+2*Ndum);
for i=1:length(f)
    kx0 = k0(i).*sin(th).*cos(ph);
    ky0 = k0(i).*sin(th).*sin(ph);    
    v = [zeros(1,Ndum) v0*exp(-1j*kx0*[0:Nx-1]*dx) zeros(1,Ndum)].';
    for n=0:Nx+2*Ndum-1
    fun = @(kx) sinc(kx*deltad/(2*pi)).^2./D_inf_Dipole_BR(kx, ky0, k0(i), dy, wd, my_lim, th, hd).*...
        exp(-1j*kx*abs(n)*dx);
    a = k0(i)/100;
    if n==0
        q = integral(fun,-50*k0(i)-1j*a,50*k0(i)+1j*a,'Waypoints',[-a-a*1j,a+a*1j]);
    else
        q = integral(fun,-a-1j*20*k0(i),k0(i)+a-1j*20*k0(i),'Waypoints',[-a-a*1j,a+a*1j, k0(i)+1j*a, k0(i)+a]);
    end
    Y1(n+1) = -q/(2*pi);
    end
    Y = toeplitz(real(Y1))+1j*toeplitz(imag(Y1));
    Z = inv(Y);
    i_delta = (Z + eye(Nx+2*Ndum)*Zl)\v;
    Znx1 = v./i_delta-Zl; 
    Znx(i,:) = Znx1;
end
figure(1);
hold on
grid on
for i=Ndum+1:length(Znx1)-Ndum
    plot(f/f0, real(Znx(:,i)),'r', 'LineWidth',2)
    plot(f/f0, imag(Znx(:,i)),'--r', 'LineWidth',2)
end
Za_dipole_BR = zeros(length(f),1);
for i=1:length(f)
    kx0 = k0(i).*sin(th).*cos(ph);
    ky0 = k0(i).*sin(th).*sin(ph);

    mx_lim = 20;
    my_lim = 20;
    Ya_dipole = zeros(size(th,1), size(th,2));
    for mx=-mx_lim:mx_lim
        kxm = kx0 - 2*pi*mx/dx;
        D_inf = D_inf_Dipole_BR(kxm, ky0, k0(i), dy, wd, my_lim, th, hd);
        Ya_dipole = Ya_dipole - (sinc(kxm*deltad/(2*pi))).^2 ./D_inf;
    end
    Ya_dipole = Ya_dipole/dx;
    Za_dipole_BR(i) = 1./Ya_dipole;
end
figure(1); plot(f/f0, real(Za_dipole_BR),'k', 'LineWidth',2)
hold on
plot(f/f0, imag(Za_dipole_BR),'--k', 'LineWidth',2)
xlabel('f/f_0')
ylabel('Re(Z_a), Im(Z_a)')
%% reflection coefficients
Gamma_f = (Znx-Zl)./(Znx+Zl);
Gamma_inf = (Za_dipole_BR-Zl)./(Za_dipole_BR+Zl);
figure()
plot(f/f0, db(abs(Gamma_inf)),'k', 'LineWidth',2)
hold on
grid on
for i=Ndum+1:length(Znx1)-Ndum
    plot(f/f0, db(abs(Gamma_f(:,i))),'r', 'LineWidth',2)
end
ylim([-25 0])
xlabel('f/f_0')
ylabel('|\Gamma|')
%%
f = 15e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;

%%

r_obs = 100000 * lambda;
l = lambda/2;
w = lambda/20;

%%

dth = pi/180; dph = dth;
% For entire region of theta and phi
[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 

% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 

%wavenumbers in all directions (x, y and z)

%kx = k0 .* sin(th) .* cos(ph);
%ky = k0 .* sin(th) .* sin(ph);
%kz = k0 .* cos(th);
ky = 0;
kx = 0:k0/50:5*k0;
kz = -1j * sqrt(-(k0^2 - kx.^2));


c = (-zeta./(2 * k0 .* kz)); %constant term in the equations

%% Calculation of the Dyad
[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

%% Calculation of Spectral Green's function

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

%% Plot of Spectral Green's function

%plot(kx, real(SGFxx))
%hold on;
%plot(kx, imag(SGFxx))

plot(kx, real(SGFzx))
hold on;
plot(kx, imag(SGFzx))
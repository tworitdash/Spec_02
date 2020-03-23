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
h = 15e-3;
%%

dth = pi/180; dph = dth;
% For entire region of theta and phi
[th, ph] = meshgrid(1e-7:dth:pi/2-dth*1e-7, pi/2); 

% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:2 * pi-dth*1e-7, pi/2); 

%wavenumbers in all directions (x, y and z)

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);


c = (-zeta./(2 * k0 * kz)); %constant term in the equations

%% Calculation of the Dyad
[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

%% Calculation of Spectral Green's function

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

%% Fourier Transform of Current distribution

%depends on the current distribution:

[Jx, Jy, Jz] = JFT_shifted_PEC(l, w, h, k0, kx, ky, kz);
[Jx_f, Jy_f, Jz_f] = JFT(l, w, k0, kx, ky, kz);

%% Far field

c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Ex_f, Ey_f, Ez_f] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx_f, Jy_f, Jz_f);

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx, Jy, Jz);

E_abs = sqrt(abs(Ex.^2 + Ey.^2 + Ez.^2));
E_abs_f = sqrt(abs(Ex_f.^2 + Ey_f.^2 + Ez_f.^2));
E_abs(91:180) = 0;

%% Directivity Calculation
U = E_abs.^2 ./ (2 * zeta);

P_rad = Prad(U, dth, dph, th, ph);

D = (4 * pi) * (U./P_rad);

%% Plots

%plot(th, D);
plot(th, E_abs/max(E_abs_f))
%hold on;
plot(th, E_abs_f/max(E_abs_f))

%surf(kx/k0, ky/k0, E_abs./max(E_abs));



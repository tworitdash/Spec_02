%%
f = 15e9;
lambda = 3e8/f;

k0 = 2 * pi / lambda;
epsilon_r = 10;
ksub = k0 * sqrt(epsilon_r);

zeta = 120*pi;

zeta_s = zeta/sqrt(epsilon_r);

%%

%%

dth = pi/180; dph = dth;
% For entire region of theta and phi
[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 

% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 

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

[Jx, Jy, Jz] = JFT(l, w, k0, kx, ky, kz);

%% Far field

c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx, Jy, Jz);

E_abs = sqrt(abs(Ex.^2 + Ey.^2 + Ez.^2));

%% Directivity Calculation
U = E_abs.^2 ./ (2 * zeta);

P_rad = Prad(U, dth, dph, th, ph);

D = (4 * pi) * (U./P_rad);

%% Plots

%plot(th, D);
%polarplot(th, E_abs/max(E_abs))

%surf(kx/k0, ky/k0, E_abs);



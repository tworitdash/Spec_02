%% Q1.2

clear all;
close all;

f = 12.5e9;
h = linspace(eps, 15e-3, 1000);
P_rad_i = zeros(size(h, 1), size(h, 2));
Dir_0 = zeros(size(h, 1), size(h, 2));

P_rad_i_d = zeros(size(h, 1), size(h, 2));
Dir_0_d = zeros(size(h, 1), size(h, 2));

for i = 1:size(h, 2)
    
lambda = 3e8./f;

zeta = 120*pi;
epsilon_r = 10;

l = 1e-3;
w = 1e-3;

r_obs = 1000 .* lambda;

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(eps:dth:pi/2-eps, eps:dph:2*pi);

k0 = 2 * pi / lambda;

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);

kz0 = -1j * sqrt(-(k0.^2 - kx.^2 - ky.^2));
k_rho_2 = kx.^2 + ky.^2;

zeta_s = zeta/sqrt(epsilon_r);

ksub = k0 * sqrt(epsilon_r);

keq = (k0 + ksub)./2;

kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));
%kzs = kz * sqrt(epsilon_r);


z = h(i)+eps;

% For TM
z0_TM = zeta * kz ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz;
zs_TE = zeta_s * ksub ./ kzs;

[zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h(i), kzs);


[vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, h(i), z, kz, kzs, z0_TE, zs_TE, z0_TM, zs_TM);


[G_xx, G_yx, G_zx] = Green(vtm, vte, itm, ite, kx, ky, k_rho_2, zeta, k0);

[Jx, Jy, Jz] = JFT_freq(l, w, keq, kx, ky);


E_far_x = 1j .* kz .* Jx .* G_xx .* exp(-1j .* kz .* abs(z - h(i))) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_y = 1j .* kz .* Jx .* G_yx .* exp(-1j .* kz .* abs(z - h(i))) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_z = 1j .* kz .* Jx .* G_zx .* exp(-1j .* kz .* abs(z - h(i))) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);


E_abs = sqrt(abs(E_far_x).^2 + abs(E_far_y).^2 + abs(E_far_z).^2);

E_th = E_far_x .* cos(th) .* cos(ph) + E_far_y .* cos(th) .* sin(ph) - E_far_z .* sin(th);
E_ph = -sin(ph) .* E_far_x + cos(ph) .* E_far_y; 


% Power radiated:

%[Dir, P_rad] = Directivity(r_obs, th, ph, E_th, E_ph, zeta, dth, dph, k0);

V_far_th = E_th * r_obs./(exp(-1j * k0 * r_obs));
V_far_ph = E_ph * r_obs./(exp(-1j * k0 * r_obs));
    
C_rad = 1/(2 * zeta);
V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);
U = C_rad * (V_far_abs).^2;

P_rad_int = C_rad * (V_far_abs).^2 .* sin(th) * dth * dph;

P_rad = sum(sum(P_rad_int));
    
Dir = 4 * pi * U ./ P_rad;

P_rad_i(i) = P_rad;
Dir_0(i) = Dir(1, 1);


% Only for one dipole in free space:
[th, ph] = meshgrid(eps:dth:pi-eps, eps:dph:2*pi);

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);

c = (-zeta./(2 * k0 .* kz)); %constant term in the equations

%% Calculation of the Dyad
[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

%% Calculation of Spectral Green's function

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

[Jx_d, Jy_d, Jz_d] = JFT_freq(l, w, k0, kx, ky);


c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx_d, Jy_d, Jz_d);

E_abs_d = sqrt(abs(Ex.^2 + Ey.^2 + Ez.^2));

E_th_d = Ex .* cos(th) .* cos(ph) + Ey .* cos(th) .* sin(ph) - Ez .* sin(th);
E_ph_d = -sin(ph) .* Ex + cos(ph) .* Ey; 


% Power radiated:

%[Dir, P_rad] = Directivity(r_obs, th, ph, E_th, E_ph, zeta, dth, dph, k0);

V_th = E_th_d * r_obs./(exp(-1j * k0 * r_obs));
V_ph = E_ph_d * r_obs./(exp(-1j * k0 * r_obs));
    
C_rad = 1/(2 * zeta);
V_abs = sqrt(abs(V_th).^2 + abs(V_ph).^2);
U_d = C_rad * (V_abs).^2;

P_rad_int_d = C_rad * (V_abs).^2 .* sin(th) * dth * dph;

P_rad_d = nansum(nansum(P_rad_int_d));
    
Dir_d = 4 * pi * U_d ./ P_rad_d;

P_rad_i_d(i) = P_rad_d;
Dir_0_d(i) = Dir_d(1, 1);

end

plot(h * 10^(3), P_rad_i./P_rad_i_d, 'LineWidth', 2);


grid on;
xlabel('height(mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('P_{rad} Normalized', 'FontSize', 12, 'FontWeight', 'bold');
title('Normalized P_{rad} with respect to a dipole radiating in free space', 'FontSize', 12, 'FontWeight', 'bold');

%print('A2_Q1_3_h', '-depsc');
%%

%%








%[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

%% Fourier Transform of Current distribution

%depends on the current distribution:



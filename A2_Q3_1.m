clear all;
close all;

f = 30e9;
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 12;
h = 5e-3;
l = lambda/2;
w = lambda/20;


r_obs = 1000 .* lambda;
dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(-pi/2-eps:dth:pi/2-eps, eps:dph:2*pi-eps);
    
k0 = 2 * pi / lambda;

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);
k_rho_2 = (kx.^2 + ky.^2);

kz0 = -1j * sqrt(-(k0.^2 - kx.^2 - ky.^2));

zeta_s = zeta/sqrt(epsilon_r);



ksub = k0 * sqrt(epsilon_r);
kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));



z = h+eps;


kxsub = k0 .* sqrt(epsilon_r) .* sin(th) .* cos(ph);
kysub = k0 .* sqrt(epsilon_r) .* sin(th) .* sin(ph);
kzsub = k0 .* sqrt(epsilon_r) .* cos(th);

k_rho_2_s = (kxsub.^2 + kysub.^2);


% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

%[zin_TE, zin_TM] = Zin(z0_TE, zs_TE, z0_TM, zs_TM, h, kz0);

[Gamma_TE, Gamma_TM] = Gamma_Q3(zs_TE, zs_TM, z0_TE, z0_TM, h, kz0);

[vte, vtm, ite, itm] = txline_3(Gamma_TE, Gamma_TM, h, kzs, kz0, z0_TE, z0_TM, zs_TE, zs_TM, z);



%Spectral Green's function

[G_xx, G_yx, G_zx] = Green_em(vtm, vte, itm, ite, kxsub, kysub, k_rho_2_s, zeta_s, ksub);

[Mx, My, Mz] = JFT_freq(l, w, ksub, kxsub, kysub);


E_far_x = 1j .* kzsub .* Mx .* G_xx .* exp(-1j .* kzsub .* abs(z)) .* exp(-1j .* ksub .* r_obs) ./ (2 * pi * r_obs);
E_far_y = 1j .* kzsub .* Mx .* G_yx .* exp(-1j .* kzsub .* abs(z)) .* exp(-1j .* ksub .* r_obs) ./ (2 * pi * r_obs);
E_far_z = 1j .* kzsub .* Mx .* G_zx .* exp(-1j .* kzsub .* abs(z)) .* exp(-1j .* ksub .* r_obs) ./ (2 * pi * r_obs);


E_abs = sqrt(abs(E_far_x).^2 + abs(E_far_y).^2 + abs(E_far_z).^2);

E_th = E_far_x .* cos(th) .* cos(ph) + E_far_y .* cos(th) .* sin(ph) - E_far_z .* sin(th);
E_ph = -sin(ph) .* E_far_x + cos(ph) .* E_far_y; 

figure(1);

plot(th(1, :)*180/pi, db(abs(E_th(91, :))/max(abs(E_th(91, :)))), 'LineWidth', 2);
grid on;
hold on;

plot(th(91, :)*180/pi, db(abs(E_ph(1, :))/max(abs(E_ph(1, :)))), 'LineWidth', 2);

ylim([-50 0]);

xlabel('\theta(deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E(\theta)_{\phi} = 90), E(\phi)_{\phi = 0})', 'FontSize', 12, 'FontWeight', 'bold');
title('Normalized Far Electric Field at different cuts', 'FontSize', 12, 'FontWeight', 'bold');
legend({['E_{\theta}(\phi = 90)'], ['E_{\phi}(\phi = 0)']},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

print('Q3_FF', '-depsc')
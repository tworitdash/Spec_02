%%
clear all;
close all;

f = 15e9;
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 10;
h = 2e-3;
l = 1e-3;
w = 1e-3;

r_obs = 1000 .* lambda;

[th, ph] = meshgrid(-pi/2-eps:pi/180:pi/2-eps, eps:pi/180:2*pi);

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


z = h+eps;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h, kzs);


[vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, h, z, kz0, kzs, z0_TE, zs_TE, z0_TM, zs_TM);


[G_xx, G_yx, G_zx] = Green(vtm, vte, itm, ite, kx, ky, k_rho_2, zeta, k0);

[Jx, Jy, Jz] = JFT_freq(l, w, keq, kx, ky);


E_far_x = 1j .* kz .* Jx .* G_xx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_y = 1j .* kz .* Jx .* G_yx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_z = 1j .* kz .* Jx .* G_zx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);


E_abs = sqrt(abs(E_far_x).^2 + abs(E_far_y).^2 + abs(E_far_z).^2);

E_th = E_far_x .* cos(th) .* cos(ph) + E_far_y .* cos(th) .* sin(ph) - E_far_z .* sin(th);
E_ph = -sin(ph) .* E_far_x + cos(ph) .* E_far_y; 

figure(1);

plot(th(1, :)*180/pi, db(E_th(1, :)/max(E_th(1, :))), 'LineWidth', 2);
hold on;
%plot(th(91, :)*180/pi, db(E_th(91, :)/max(E_th(91, :))), 'LineWidth', 2);
grid on;

hold on;

plot(th(91, :)*180/pi, db(E_ph(91, :)/max(E_ph(91, :))), 'LineWidth', 2);
hold on;
%plot(th(1, :)*180/pi, db(E_ph(1, :)/max(E_ph(1, :))), 'LineWidth', 2);

xlabel('\theta(deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E(\theta)_{\phi} = 0), E(\phi)_{\phi = 90})', 'FontSize', 12, 'FontWeight', 'bold');
title('Normalized Far Electric Field at different cuts', 'FontSize', 12, 'FontWeight', 'bold');
legend({'E_{\theta}(\phi = 0)', 'E_{\phi}(\phi = 90)'}, 'Location', 'south', 'FontSize', 12, 'FontWeight', 'bold');

ylim([-40 0]);
print('Q1_Far_field', '-depsc');
figure(2);

plot(th(1, :)*180/pi, db(E_abs(1, :)/max(E_abs(1, :))), '--', 'LineWidth', 2);
hold on;
plot(th(91, :)*180/pi, db(E_abs(91, :)/max(E_abs(91, :))), '*', 'LineWidth', 1);
grid on;


xlabel('\theta(deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E(\theta, \phi)', 'FontSize', 12, 'FontWeight', 'bold');
title('Far field', 'FontSize', 12, 'FontWeight', 'bold');

grid on;
ylim([-40 0])

%print('Q1_Spec_Gej', '-depsc');



%%

%%








%[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

%% Fourier Transform of Current distribution

%depends on the current distribution:



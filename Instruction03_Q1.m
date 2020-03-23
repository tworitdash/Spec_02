%% 
clear all;
close all;

f = 20e9;
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 2.2;
h = (3e8/10e9/sqrt(epsilon_r))/4;
l = (3e8/10e9)/2;
w = 1e-3;

r_obs = 1000 .* lambda;

[th, ph] = meshgrid(-pi/2+eps:pi/180:pi/2-eps, eps:pi/180:2*pi);

k0 = 2 * pi ./ lambda;
ksub = k0 * sqrt(epsilon_r);

% kx = linspace(k0, ksub, 100);
% ky = 0;


%k_rho_2 = kx.^2 + ky.^2;

krho_i = linspace(k0, ksub, 100);

kz0 = -1j * sqrt(-(k0.^2 - krho_i.^2));

zeta_s = zeta/sqrt(epsilon_r);

keq = (k0 + ksub)./2;

kzs = -1j * sqrt(-(ksub.^2 - krho_i.^2));


z = h+eps;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h, kzs);

[D_krho_TE] = D_kro(k0, "TE", epsilon_r, h, zeta, zeta_s, krho_i);
[D_krho_TM] = D_kro(k0, "TM", epsilon_r, h, zeta, zeta_s, krho_i);

% D_krho_TE = zup_TE + zdn_TE;
% D_krho_TM = zup_TM + zdn_TM;

figure(1);

plot(krho_i/k0, abs(1./D_krho_TE), 'LineWidth', 2);
hold on;
plot(krho_i/k0, abs(1./D_krho_TM), 'LineWidth', 2);
grid on;

xlabel('k_{\rho}/k_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('1./Denominator D(k_{\rho})', 'FontSize', 12, 'FontWeight', 'bold');
title(['Denominator function for finding poles at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'D(TE)', 'D(TM)'}, 'Location', 'north', 'FontSize', 12, 'FontWeight', 'bold');

print(['A3_Q1_K_guess_' , num2str(epsilon_r)], '-depsc');

%krho_g_TE = 1.874 .* k0;
%krho_g_TM = 2.616 .* k0;
% krho_g_TE = 1.025 .* k0;
% krho_g_TM = 1.587 .* k0;

krho_g_TE = 1.122 .* k0;
krho_g_TM = 1.356 .* k0;

% krho_g_TE = 0;
% krho_g_TM = 1.112 .* k0;

[krho_TE, f_axis_1] = finddrop(k0, epsilon_r, krho_g_TE, h, zeta, zeta_s, "TE");
[krho_TM, f_axis_2] = finddrop(k0, epsilon_r, krho_g_TM, h, zeta, zeta_s, "TM");



figure(2);

plot(f_axis_1 * 10^(-9), abs(krho_TE), 'LineWidth', 2);
hold on;
plot(f_axis_2 * 10^(-9), abs(krho_TM), 'LineWidth', 2);
grid on;

%ylim([1 2]);

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('k_{\rho g}(f_i)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Guess pole point vs Frequency  at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'k_{\rho g}(f_i)(TE)', 'k_{\rho g}(f_i)(TM)'}, 'Location', 'north', 'FontSize', 12, 'FontWeight', 'bold');

print(['A3_Q1_K_guess_2_', num2str(epsilon_r)], '-depsc');

%% 

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



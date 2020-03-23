%% 
clear all;
close all;

f = 20e9;
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 10;
h = 2e-3;


r_obs = 1000 .* lambda;

[th, ph] = meshgrid(-pi/2-eps:pi/180:pi/2-eps, eps:pi/180:2*pi);

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


z = h/2;

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

plot(krho_i/k0, real(D_krho_TE), 'LineWidth', 2);
hold on;
plot(krho_i/k0, imag(D_krho_TE), 'LineWidth', 2);
hold on;
plot(krho_i/k0, abs(D_krho_TE), 'LineWidth', 2);
hold on;
% plot(krho_i/k0, abs(D_krho_TM), 'LineWidth', 2);
% grid on;

xlabel('k_{\rho}/k_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Denominator D(k_{\rho})', 'FontSize', 12, 'FontWeight', 'bold');
title(['Denominator function for finding poles at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'D(TE)', 'D(TM)'}, 'Location', 'north', 'FontSize', 12, 'FontWeight', 'bold');

%print(['A3_Q1_K_guess_' , num2str(epsilon_r)], '-depsc');

krho_g_TE = 1.874 .* k0;
krho_g_TM = 2.616 .* k0;

[krho_TE, f_axis_1] = finddrop(k0, epsilon_r, krho_g_TE, h, zeta, zeta_s, "TE");
[krho_TM, f_axis_2] = finddrop(k0, epsilon_r, krho_g_TM, h, zeta, zeta_s, "TM");



figure(2);

plot(f_axis_1 * 10^(-9), abs(krho_TE), 'LineWidth', 2);
hold on;
plot(f_axis_2 * 10^(-9), abs(krho_TM), 'LineWidth', 2);
grid on;

ylim([1 3]);

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('k_{\rho g}(f_i)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Guess pole point vs Frequency  at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'k_{\rho g}(f_i)(TE)', 'k_{\rho g}(f_i)(TM)'}, 'Location', 'north', 'FontSize', 12, 'FontWeight', 'bold');

%print(['A3_Q1_K_guess_2_', num2str(epsilon_r)], '-depsc');

f0 = 10e9;
lambda0 = 3e8./f0;

k0_ = (2*pi*f0)./(3e8);

k_sw = krho_TM(f_axis_1==f0) .* k0_;



[VtmR, ItmR] = Residue_GroundSlab(k0_, epsilon_r, h, k_sw, z);

[rho, phi] = meshgrid(1:0.01:100, eps:pi/180:2*pi);



[Erho, Ez, Hphi] = TMSwFields(k0_, k_sw, epsilon_r, VtmR, ItmR, rho, phi, z, h);

figure(3);

% plot(rho(1, :), real(Erho(1, :)), 'LineWidth', 2);
% hold on;
% plot(rho(1, :), imag(Erho(1, :)), '--', 'LineWidth', 2);
% hold on;
% plot(rho(1, :), real(Ez(1, :)), 'LineWidth', 2);
% hold on;
% plot(rho(1, :), imag(Ez(1, :)), '--', 'LineWidth', 2);
% 
% %xlim([1 ]);
% 
% xlabel('\rho [m]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('E_{\rho}, E_{z}', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Far field of surface waves inside the substrate'], 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Real( E_{\rho} )', 'Imag( E_{\rho} )', 'Real( E_{z} )', 'Imag( E_{z})'}, 'Location', 'north', 'FontSize', 12, 'FontWeight', 'bold');


plot(phi(:, 1)*180/pi, abs(Erho(:, 1)), 'LineWidth', 2, 'color', [0.6350 0.0780 0.1840]);

hold on;
plot(phi(:, 1)*180/pi, abs(Ez(:, 1)), 'LineWidth', 2, 'color', [0.4940 0.1840 0.5560]);


grid on;

xlabel('\phi(deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_{\rho}, E_{z}', 'FontSize', 12, 'FontWeight', 'bold');
title(['Far field of surface waves inside the substrate vs \phi'], 'FontSize', 12, 'FontWeight', 'bold');
legend({'E_{\rho}', 'E_{z}'}, 'Location', 'southeast', 'FontSize', 12, 'FontWeight', 'bold');

%print('A3_Q2_FF_phi', '-depsc');

figure(4);

polarplot(phi(:, 1), abs(Erho(:, 1)), 'LineWidth', 2, 'color', [0.6350 0.0780 0.1840]);

hold on;
polarplot(phi(:, 1), abs(Ez(:, 1)), 'LineWidth', 2, 'color', [0.4940 0.1840 0.5560]);
legend({'E_{\rho}', 'E_{z}'}, 'Location', 'southeast', 'FontSize', 12, 'FontWeight', 'bold');

grid on;




%print('A3_Q2_FF_phi_polar', '-depsc');






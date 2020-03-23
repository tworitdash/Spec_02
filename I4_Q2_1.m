%%
close all;
clear all;

f = 31e9;
c = 3e8;
lambda = c./f;
zeta = 120*pi;
epsilon_r = linspace(1, 25, 50);
h = 5e-3;


krho_TE_e = zeros(size(epsilon_r));
krho_TM_e = zeros(size(epsilon_r));

krho_TE_a = zeros(size(epsilon_r));
krho_TM_a = zeros(size(epsilon_r));
 
for m = 1:size(epsilon_r, 2)

l = lambda./2;
w = lambda./20;

k0 = 2 * pi / lambda;

f1 = 30e9;
lambda1 = c./f1;
k1 = (2 .* pi)./lambda1;

% hs = lambda1./(4 .* sqrt(epsilon_r(m))); No hs for resonant lens feed antenna

zeta_s = zeta/sqrt(epsilon_r(m));

ksub = k0 * sqrt(epsilon_r(m));

[krho_real, krho_imag] = meshgrid(linspace(eps, k0, 1000), linspace(eps, k0, 1000));

krho = krho_real - 1j .* krho_imag;

[D_krho_TE] = D_kro_leak_lens(k0, "TE", epsilon_r(m), h, zeta, zeta_s, krho);
[D_krho_TM] = D_kro_leak_lens(k0, "TM", epsilon_r(m), h, zeta, zeta_s, krho);

D_TE_max = max(abs(1./D_krho_TE(:)));
D_TM_max = max(abs(1./D_krho_TM(:)));

[Ind_x_TE, Ind_y_TE] = find(abs(1./D_krho_TE) == D_TE_max);
[Ind_x_TM, Ind_y_TM] = find(abs(1./D_krho_TM) == D_TM_max);

krho_g_TE = krho(Ind_x_TE, Ind_y_TE);

krho_g_TM = krho(Ind_x_TM, Ind_y_TM);

% figure(9);
% contourf(krho_real/k0, -krho_imag/k0, abs(1./D_krho_TE));
% figure(10);
% contourf(krho_real/k0, -krho_imag/k0, abs(1./D_krho_TM));

%% finddrop

[krho_TE, f_axis_1] = finddrop_I4_lens(k0, epsilon_r(m), krho_g_TE, h, zeta, zeta_s, "TE");
[krho_TM, f_axis_2] = finddrop_I4_lens(k0, epsilon_r(m), krho_g_TM, h, zeta, zeta_s, "TM");

krho_TE_e(m) = krho_TE(find(abs(f_axis_1 - f1) < 1e-5)) .* k1;
krho_TM_e(m) = krho_TM(find(abs(f_axis_1 - f1) < 1e-5)) .* k1;

%% Analytic solutions of Dispertion equations:



hb = h./lambda1;


kz0_TE = (k1./(pi.*sqrt(epsilon_r(m)).*(2.*hb).^2)) .* ((2.*pi.*hb.*sqrt(epsilon_r(m)) + 1j)./(1 + (1./(2.*pi.*hb.*epsilon_r(m)).^2)));
kz0_TM = (k1./(4 .* hb)).*(1 + sqrt(1 + 8j .* (hb./(pi.*sqrt(epsilon_r(m))))));

krho_TE_a(m) = sqrt(k1.^2 - kz0_TE.^2);

krho_TM_a(m) = sqrt(k1.^2 - kz0_TM.^2);






end
figure(1);

plot(epsilon_r, real(krho_TE_e./k1), 'LineWidth', 2);
hold on;
plot(epsilon_r, real(krho_TE_a./k1), '--', 'LineWidth', 2);
hold on;
plot(epsilon_r, real(krho_TM_e./k1), 'LineWidth', 2);
hold on;
plot(epsilon_r, real(krho_TM_a./k1), '--', 'LineWidth', 2);

grid on;

xlim([1.49 25]);

xlabel('\epsilon_r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('k_{\rho g}(f_i)/k_0', 'FontSize', 12, 'FontWeight', 'bold');
title(['Guess pole point vs \epsilon_r'], 'FontSize', 12, 'FontWeight', 'bold');
legend({'Re(k_{\rho g}(f_i)(TE))', 'Re(k_{\rho g}(f_i)(TE)) Analytical', 'Re(k_{\rho g}(f_i)(TM))', 'Re(k_{\rho g}(f_i)(TM)) Analytical'}, 'Location', 'northeast', 'FontSize', 12, 'FontWeight', 'bold');
print('A4_Q2_real_er', '-depsc');




figure(2);

plot(epsilon_r, imag(krho_TE_e./k1), 'LineWidth', 2);
hold on;
plot(epsilon_r, imag(krho_TE_a./k1), '--', 'LineWidth', 2);
hold on;
plot(epsilon_r, imag(krho_TM_e./k1), 'LineWidth', 2);
hold on;
plot(epsilon_r, imag(krho_TM_a./k1), '--', 'LineWidth', 2);

xlim([1.49 25]);

grid on;

xlabel('\epsilon_r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('k_{\rho g}(f_i)/k_0', 'FontSize', 12, 'FontWeight', 'bold');
title(['Guess pole point vs \epsilon_r'], 'FontSize', 12, 'FontWeight', 'bold');
legend({'Im(k_{\rho g}(f_i)(TE))', 'Im(k_{\rho g}(f_i)(TE)) Analytical', 'Im(k_{\rho g}(f_i)(TM))', 'Im(k_{\rho g}(f_i)(TM)) Analytical'}, 'Location', 'southeast', 'FontSize', 12, 'FontWeight', 'bold');
print('A4_Q2_imag_er', '-depsc');

%z = h + hs;

% [zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM] = Zall(z0_TE, zs_TE, z0_TM, zs_TM, h, hs, kzs, kz0);
% 
% [Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM] = Gamma(zs_TE, zs_TM, zL_TE, zL_TM, z0_TE, z0_TM, h, hs, kz0, kzs);
% 
% [vte, vtm, ite, itm] = txline_2(Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM, h, hs, kzs, kz0, zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM, z);